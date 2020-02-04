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



    #include "mars_constant.hpp"
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
					std::cout<<it<<"/"<<max_iter<<std::endl;

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
		class DenseMatrix
		{
		public:
			using Vec=std::vector<Real>;
			using VecVec=std::vector<std::vector<Real>>;
			


			// inline void init(std::vector<std::vector<Real>>& mat)
			// {
			// 	// A_ptr_->clear();
			// 	A_ptr_=std::make_shared<VecVec>(mat);
			// 	rows_=mat.size();
			// 	// we assume all the rows have the same number of columns
			// 	cols_=mat[0].size();

			// }
            DenseMatrix()        
            {}

            DenseMatrix(const Integer rows, const Integer cols):
 			rows_(rows),
			cols_(cols)          
            {vec_.resize(rows*cols);}

			inline void init(const Integer rows, const Integer cols)
			{
			  vec_.resize(rows*cols);
			  rows_=rows;
			  cols_=cols;
			}

			inline void resize(const Integer rows, const Integer cols)
			{
			  vec_.resize(rows*cols);
			  rows_=rows;
			  cols_=cols;
			}

			inline auto& get(const Integer i, const Integer j)
			{
			  return vec_[i*cols_+j];
			}

			inline auto& get(const Integer i, const Integer j)const 
			{
			  return vec_[i*cols_+j];
			}

			inline auto& operator()(const Integer i, const Integer j)
			{
			  // std::cout<<"(i,j)"<<i<<","<<j<<std::endl;
			  return vec_[i*cols_+j];
			}

			inline auto& operator()(const Integer i, const Integer j)const 
			{
			  // std::cout<<"(i,j)"<<i<<","<<j<<std::endl;
			  return vec_[i*cols_+j];
			}


			// inline void cholesky_matrix(VecVec& Chol)
			// {
			// 	if(Chol.size()==0)
			// 		Chol.resize(rows_);
			// 	const auto& A=*A_ptr_;
			// 	Chol[0].resize(1);
   //              Chol[0][0]=sqrt(A[0][0]);
			// 	for(Integer i=1;i<rows_;i++)
			// 		{
			// 			Chol[i].resize(i+1);
			// 			for(Integer j=0;j<i;j++)
			// 			{
			// 			  Chol[i][j]=A[i][j];
			// 			  for(Integer k=0;k<j;k++)
			// 			  	Chol[i][j]-=Chol[i][k]*Chol[j][k];

			// 			  Chol[i][j]=Chol[i][j]/Chol[j][j];
			// 			}

			// 			Chol[i][i]=A[i][i];
			// 			for(Integer k=0;k<i;k++)
			// 		       Chol[i][i]-=Chol[i][k]*Chol[i][k];
   //                       Chol[i][i]=sqrt(Chol[i][i]);
			// 		}

			// }

			inline void cholesky_matrix(DenseMatrix& Chol)
			{
       
                Chol(0,0)=sqrt(get(0,0));



				for(Integer i=1;i<rows_;i++)
					{
						// std::cout<<" "<<std::endl;
						// std::cout<<"row i="<<i<<std::endl;
						for(Integer j=0;j<i;j++)
						{
							// std::cout<<" "<<std::endl;
							// std::cout<<"col j="<<j<<std::endl;
						  Chol(i,j)=get(i,j);
						  // std::cout<<Chol(i,j)<<" ";
						  for(Integer k=0;k<j;k++)
						  	{
						  		Chol(i,j)-=Chol(i,k)*Chol(j,k);
						  		// std::cout<<Chol(i,j)<<" ";
						  	}
						  // std::cout<<" "<<std::endl;
                          // std::cout<<"Chol("<<j<<","<<j<<")="<<Chol(j,j)<<std::endl;
						  Chol(i,j)=Chol(i,j)/Chol(j,j);
						  // std::cout<<" "<<std::endl;
						  // std::cout<<"final= "<<std::endl;
						  // std::cout<<Chol(i,j)<<" ";
						  // std::cout<<" "<<std::endl;


		                // std::cout<<"MAT "<<std::endl;
		                // for(Integer s=0;s<rows_;s++)
		                // {

		                // 	for(Integer t=0;t<cols_;t++)
		                // 		std::cout<<Chol(s,t)<<" ";
		                // 	std::cout<<" "<<std::endl;
		                // }
		                // std::cout<<" "<<std::endl;


						}

						Chol(i,i)=get(i,i);
						// std::cout<<" "<<std::endl;
						// std::cout<<Chol(i,i)<<" ";
						for(Integer k=0;k<i;k++)
					       {
					       	Chol(i,i)-=Chol(i,k)*Chol(i,k);
					       	// std::cout<<Chol(i,i)<<" ";
					       }
                         Chol(i,i)=sqrt(Chol(i,i));

                         // std::cout<<" "<<std::endl;

                         // std::cout<<"col j="<<i<<std::endl;

                         // std::cout<<Chol(i,i)<<" ";
                         // std::cout<<std::endl;

						for(Integer j=i+1;j<cols_;j++)
							Chol(i,j)=0;
		                // std::cout<<"MAT "<<std::endl;
		                // for(Integer s=0;s<rows_;s++)
		                // {

		                // 	for(Integer t=0;t<cols_;t++)
		                // 		std::cout<<Chol(s,t)<<" ";
		                // 	std::cout<<" "<<std::endl;
		                // }
		                // std::cout<<" "<<std::endl;
						

					}

			}


			inline void transpose_upper_triangular_solve(Vec& x,const VecVec& L,const Vec& b)
			{
				x[0]=b[0]/L[0][0];

				for(Integer i=1;i<rows_;i++)
				{   
					x[i]=b[i];

					for(Integer k=0;k<i;k++)
						x[i]-=L[k][i]*x[k];

					x[i]=x[i]/L[i][i];
				}
			}

			// inline void transpose_lower_triangular_solve(Vec& x,const VecVec& L,const Vec& b)
			// {
   //              Integer n=rows_-1;
			// 	x[n]=b[n]/L[n][n];

			// 	for(Integer i=n-1;i>=0;i--)
			// 	{   

			// 		x[i]=b[i];

			// 		for(Integer k=i+1;k<n+1;k++)
			// 			x[i]-=L[k][i]*x[k];

			// 		x[i]=x[i]/L[i][i];
			// 		std::cout<<"x[i]="<<x[i]<<std::endl;
			// 	}

			// }
			// inline void lower_triangular_solve(Vec& x,const VecVec& L,const Vec& b)
			// {
			// 	x[0]=b[0]/L[0][0];

			// 	for(Integer i=1;i<rows_;i++)
			// 	{   
			// 		x[i]=b[i];

			// 		for(Integer k=0;k<i;k++)
			// 			x[i]-=L[i][k]*x[k];

			// 		x[i]=x[i]/L[i][i];
			// 	}
			// }
 
 			inline void transpose_lower_triangular_solve(Vec& x,const DenseMatrix& L,const Vec& b)
			{
                Integer n=rows_-1;
				x[n]=b[n]/L(n,n);

				for(Integer i=n-1;i>=0;i--)
				{   

					x[i]=b[i];

					for(Integer k=i+1;k<n+1;k++)
						x[i]-=L(k,i)*x[k];

					x[i]=x[i]/L(i,i);
					// std::cout<<"x[i]="<<x[i]<<std::endl;
				}

			}
			inline void lower_triangular_solve(Vec& x,const DenseMatrix& L,const Vec& b)
			{
				x[0]=b[0]/L(0,0);

				for(Integer i=1;i<rows_;i++)
				{   
					x[i]=b[i];

					for(Integer k=0;k<i;k++)
						x[i]-=L(i,k)*x[k];

					x[i]=x[i]/L(i,i);
				}
			}

 			inline void upper_triangular_solve(Vec& x,const VecVec& L,const Vec& b)
			{
                Integer n=rows_-1;
				x[n]=b[n]/L[n][n];

				for(Integer i=n-1;i>=0;i--)
				{   
					x[i]=b[i];

					for(Integer k=i+1;k<n;k++)
						x[i]-=L[i][k]*x[k];

					x[i]=x[i]/L[i][i];
				}
			}
			inline void cholesky(std::vector<Real>&x,const std::vector<Real>&b)
			{
				// we solve Ax=b as H H' x = b, so H'y=b and H x=y

				// std::cout<< "pre cholesky"<<std::endl;
				Vec y(rows_);
				DenseMatrix<Real> H;
				// std::cout<<"rows cols=="<<rows_<<", "<<cols_<<std::endl;
				H.resize(rows_,cols_);
				cholesky_matrix(H);
    //             std::cout<<" "<<std::endl;
				// std::cout<<"chol"<<std::endl;

				// for(Integer i=0;i<rows_;i++)
				// {
				// 	for(Integer j=0;j<=cols_;j++)
				// 	std::cout<<H(i,j)<<" ";
				//  std::cout<<std::endl;
				// }

				lower_triangular_solve(y,H,b);

				// std::cout<<"y  "<<std::endl;

				// for(Integer i=0;i<y.size();i++)
				// 	std::cout<<y[i]<<std::endl;
				// std::cout<<std::endl;

				// std::cout<<"x  "<<std::endl;

				transpose_lower_triangular_solve(x,H,y);
				// for(Integer i=0;i<x.size();i++)
				// 	std::cout<<x[i]<<std::endl;
				// std::cout<<std::endl;

				// std::cout<< "after cholesky"<<std::endl;
				// upper_triangular_solve(x,H,y);


			}

			Integer rows(){return rows_;}
			Integer cols(){return cols_;}



		private:
			Vec vec_;
			Integer rows_;
			Integer cols_;
		};





	template<typename T>
		class SparseMatrix
		{
		private:
			std::vector<T> A_;
			std::vector<std::map<Integer,Integer>> cols_idx_; 
			std::vector<std::map<Integer,Integer>> rows_idx_; 
			Integer max_rows_;
			Integer max_cols_;
			Integer max_cols_used_;
			Integer non_zeros_; 


		public:
			SparseMatrix(){}


			SparseMatrix(const Integer max_rows,const Integer max_cols):
			max_rows_(max_rows),
			max_cols_(max_cols),
			non_zeros_(0)
			{
				A_.resize(max_rows_*max_cols);
				cols_idx_.resize(max_rows_);
				rows_idx_.resize(max_cols_);
			}

			auto rows(){return max_rows_;}
			void init(const Integer max_rows,const Integer max_cols,const Integer max_cols_used)
			{
				max_rows_=max_rows;
				max_cols_=max_cols;
				max_cols_used_=max_cols_used;
				// std::cout<<"SparseMatrix"<<std::endl;
				// std::cout<<"max_row="<<max_rows<<std::endl;
				// std::cout<<"max_cols="<<max_cols<<std::endl;
				// std::cout<<"max_cols_used="<<max_cols_used<<std::endl;

				non_zeros_=0;
				A_.resize(max_rows_*max_cols_used);
				cols_idx_.resize(max_rows_);
				rows_idx_.resize(max_cols_);
			}

			void equal(const T& value,const Integer i, const Integer j)
			{
				// std::cout<<"(i,j)==("<<i<<","<<j<<")"<<std::endl;
				if(cols_idx_[i].count(j))
				{
                    
                    // std::cout<<"already exists in "<< cols_idx_[i].count(j)<<std::endl;
                    // std::cout<<"previous  "<< A_[cols_idx_[i].at(j)]<<std::endl;
                    // std::cout<<"new  "<< value<<std::endl;

					A_[cols_idx_[i].at(j)]=value;

					
					// for(std::size_t i=0;i<A_.size();i++)
					//   std::cout<<A_[i]<<" ";
					// std::cout<<std::endl;
				}
				else
				{
					
					// std::cout<<"does not already exists"<<std::endl;
					A_[non_zeros_]=value ;
					// for(std::size_t i=0;i<A_.size();i++)
					// std::cout<<A_[i]<<" ";
					// std::cout<<std::endl;
					cols_idx_[i].insert (std::pair<Integer,Integer>(j,non_zeros_) );
					rows_idx_[j].insert (std::pair<Integer,Integer>(i,non_zeros_) );
					non_zeros_++ ; 
                    

     //                std::cout<<"["<<i<<","<<j<<"]"<<std::endl;
     //                std::cout<<"cols_idx_"<<std::endl;
					// for (auto it=cols_idx_[i].begin(); it!=cols_idx_[i].end(); ++it)
					// 	std::cout<<"( " << i<<", "<< it->first << ", " << A_[it->second]<<" )" ;
					// std::cout<<std::endl;

					// std::cout<<"rows_idx_"<<std::endl;
					// for (auto it=rows_idx_[j].begin(); it!=rows_idx_[j].end(); ++it)
					// 	std::cout<<"( " << i<<", "<< it->first << ", " << A_[it->second]<<" )" ;
					// std::cout<<std::endl;


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
					A_[non_zeros_]=value ;

// for(std::size_t i=0;i<A_.size();i++)
//   std::cout<<A_[i]<<" ";
// std::cout<<std::endl;
					cols_idx_[i].insert (std::pair<Integer,Integer>(j,non_zeros_) );
					rows_idx_[j].insert (std::pair<Integer,Integer>(i,non_zeros_) );
					non_zeros_++ ; 
				}

			}




			void set_zero_row(const Integer i)
			{
				const auto& map=cols_idx_[i]; 
				for (auto it=map.begin(); it!=map.end(); ++it)
					A_[it->second]=0;
			}

			void print_row(const Integer i)
			{
				// std::cout<<"printing matrix row "<<i<<std::endl;

					const auto& map=cols_idx_[i]; 
					for (auto it=map.begin(); it!=map.end(); ++it)
						std::cout<<"( " << i<<", "<< it->first << ", " << A_[it->second]<<" )" ;
					std::cout<<std::endl;
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

			void print_row_val(const Integer i)const
			{
				Integer cont_cols=0;
				// std::cout<<"printing matrix row "<<i<<std::endl;
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
				for (int s=cont_cols;s<max_cols_;s++)
					std::cout<<0<<" ";

					std::cout<<std::endl;
			}

			void print_val()const
			{
				// Integer cont_cols=0;

				std::cout<<"printing matrix"<<std::endl;
				for(Integer i=0;i<max_rows_;i++)
				{   
					print_row_val(i);
				}
			}

			auto multiply(const Integer i,const std::vector<T>& b)const
			{
				T tmp=0;
				const auto& map=cols_idx_[i]; 
// std::cout<<" multiply ";
				// std::cout<<"row="<<i<<std::endl;
				for (auto it=map.begin(); it!=map.end(); ++it)
				{
					// std::cout<<"j=="<<it->first<<" A(i,j)="<<A_[it->second]<<" b(j)="<< b[it->first]<< std::endl;
					tmp+=A_[it->second]*b[it->first];
// std::cout<<it->first<<"  "<<it->second<<"  A="<<A_[it->second]<<"   b="<<b[it->first]<<" "<<std::endl;
// std::cout<<tmp<<" "<<std::endl;;
				}
                 // std::cout<<tmp<<std::endl;
				return tmp;
			}

			void multiply(std::vector<T>& result,const std::vector<T>& b)const
			{
				// std::cout<< " multiply result=A*b"<<std::endl;
				if(result.size()==0)
					result.resize(max_rows_);
				for(Integer i=0;i<max_rows_;i++)
				{
					result[i]=multiply(i,b);
				}
			}

			void multiply(std::vector<T>& result,const std::vector<T>& b,const std::vector<Integer>& index)const
			{
				// std::cout<< " multiply result=A*b"<<std::endl;
				// if(result.size()==0)
					result.resize(index.size());
				for(Integer i=0;i<index.size();i++)
				{
					result[i]=multiply(index[i],b);
				}
			}

			auto transpose_multiply(const Integer i,const std::vector<T>& b) const
			{
				T tmp=0;
				const auto& map=rows_idx_[i]; 
// std::cout<<" multiply ";
				// std::cout<<"row="<<i<<std::endl;
				for (auto it=map.begin(); it!=map.end(); ++it)
				{
					// std::cout<<"j=="<<it->first<<" A(i,j)="<<A_[it->second]<<" b(j)="<< b[it->first]<< std::endl;
					tmp+=A_[it->second]*b[it->first];
// std::cout<<it->first<<"  "<<it->second<<"  A="<<A_[it->second]<<"   b="<<b[it->first]<<" "<<std::endl;
// std::cout<<tmp<<" "<<std::endl;;
				}
                 // std::cout<<tmp<<std::endl;
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

			void row_static_condensation(const Integer i,const std::vector<bool>& constraint)
			{
				const auto& map=cols_idx_[i]; 
// std::cout<<" multiply ";
				for (auto it=map.begin(); it!=map.end(); ++it)
				{
					if(constraint[it->first])
					{
						A_[it->second]=0;
					}
// std::cout<<it->first<<"  "<<it->second<<"  A="<<A_[it->second]<<"   b="<<b[it->first]<<" "<<std::endl;
// std::cout<<tmp<<" "<<std::endl;;
				}
			}

			auto multiply(std::vector<T>& b)
			{

				std::vector<T> result(max_rows_);
				// std::cout<<"multiply="<<std::endl;
// std::cout<<std::endl<<std::endl;
				for(Integer i=0;i<max_rows_;i++)
				{
					result[i]=multiply(i,b);
// std::cout<<std::endl<<result[i]<<std::endl;
				}
// std::cout<<std::endl<<std::endl;
				return result;
			}

			void multiply_and_add(std::vector<T>& result,const Real alpha,const std::vector<T>& x,const std::vector<T>& b, const Real beta=1.0 )const
			{
                if(result.size()==0)
				result.resize(max_rows_);
			// std::cout<<"max rows=="<<max_rows_<<std::endl;
			// std::cout<<"b.size()=="<<b.size()<<std::endl;
			// std::cout<<"x.size()=="<<x.size()<<std::endl;
				for(Integer i=0;i<max_rows_;i++)
				{
					result[i]=alpha * multiply(i,x)+ b[i]*beta;
				}
				// return result;
			}




			auto multiply(const SparseMatrix& B)
			{

				SparseMatrix C;
				// std::cout<<"mat mat multiply" <<std::endl;

			    C.set_max_rows(max_rows_);

			    C.set_max_cols(B.max_cols());

			    auto& B_cols_idx=B.cols_idx();

			    auto& B_rows_idx=B.rows_idx();

			    auto B_max_cols=B.max_cols();

				C.cols_idx().resize(C.max_rows());
				C.rows_idx().resize(C.max_cols());


			    std::map<Integer,Integer> intersection;
                Real tmp;
                C.non_zeros()=0;
 
				for(Integer i=0;i<max_rows_;i++)
				{
					for(Integer j=0;j<B_max_cols;j++)
					{


						if(!intersection.empty())
						{
						intersection.clear();
						}						
                    std::set_intersection(cols_idx_[i].begin(),cols_idx_[i].end(),
                    					  B_rows_idx[j].begin(),B_rows_idx[j].end(),
									      inserter(intersection, intersection.begin()),
									      [](const std::pair<Integer,Integer>& p1,const std::pair<Integer,Integer>& p2){
									         return p1.first < p2.first;});
                    if(!intersection.empty())
                    {
                    	tmp=0.0;                    	
						for (auto it=intersection.begin(); it!=intersection.end(); ++it)
						{ 
							tmp+=operator()(i,it->first)*B(it->first,j);
						}
                        C().push_back(tmp);
						C.cols_idx()[i].insert (std::pair<Integer,Integer>(j,C.non_zeros()) );
						C.rows_idx()[j].insert (std::pair<Integer,Integer>(i,C.non_zeros()) );
						C.non_zeros()++ ; 
                    }
					}
				}
				return C;
			}


			auto transpose_and_multiply(std::vector<T>& b)
			{
				std::vector<T> result(max_cols_);
				for(Integer i=0;i<max_cols_;i++)
				{
					result[i]=transpose_multiply(i,b);
				}
				return result;
			}

			auto transpose_and_multiply(std::vector<T>& result,const std::vector<T>& b)const
			{
				if(result.size()!=max_cols_)
					result.resize(max_cols_);
				for(Integer i=0;i<max_cols_;i++)
				{
					result[i]=transpose_multiply(i,b);
				}
				return result;
			}

			auto transpose_and_multiply(const SparseMatrix& B)
			{

				SparseMatrix C;
				std::cout<<"mat transpose and mat multiply" <<std::endl;

			    C.set_max_rows(max_cols_);

			    C.set_max_cols(B.max_cols());

			    auto& B_cols_idx=B.cols_idx();

			    auto& B_rows_idx=B.rows_idx();

			    auto B_max_cols=B.max_cols();

				C.cols_idx().resize(C.max_rows());
				C.rows_idx().resize(C.max_cols());


			    std::map<Integer,Integer> intersection;
                Real tmp;
                C.non_zeros()=0;
 
				for(Integer i=0;i<max_cols_;i++)
				{
					for(Integer j=0;j<B_max_cols;j++)
					{


					if(!intersection.empty())
					{
					intersection.clear();
					}						
                    std::set_intersection(rows_idx_[i].begin(),rows_idx_[i].end(),
                    					  B_rows_idx[j].begin(),B_rows_idx[j].end(),
									      inserter(intersection, intersection.begin()),
									      [](const std::pair<Integer,Integer>& p1,const std::pair<Integer,Integer>& p2){
									         return p1.first < p2.first;});

                    // std::cout<< "before entering in (i,j)=="<<i<<","<<j<<std::endl;
                    if(!intersection.empty())
                    {
                    	// std::cout<< "after entering in (i,j)=="<<i<<","<<j<<std::endl;
                    	tmp=0.0;                    	
						for (auto it=intersection.begin(); it!=intersection.end(); ++it)
						{ 
							tmp+=operator()(it->first,i)*B(it->first,j);
						}
                        C().push_back(tmp);
						C.cols_idx()[i].insert (std::pair<Integer,Integer>(j,C.non_zeros()) );
						C.rows_idx()[j].insert (std::pair<Integer,Integer>(i,C.non_zeros()) );
						C.non_zeros()++ ; 
                    }
					}
				}
				return C;
			}




            inline auto& operator()() {return A_;}

			auto&  operator()(const Integer i, const Integer j)
			{
				return A_[cols_idx_[i].at(j)];
			}


			auto&  operator() (const Integer i, const Integer j)const
			{
				return A_[cols_idx_[i].at(j)];
			}

			inline auto&  get(const Integer i, const Integer j)
			{
				return A_[cols_idx_[i].at(j)];
			}

			inline auto&  get(const Integer i, const Integer j)const
			{
				return A_[cols_idx_[i].at(j)];
			}

			inline Real  get_element_or_zero(const Integer i, const Integer j)const
			{
				if(cols_idx_[i].count(j))
					return A_[cols_idx_[i].at(j)];
			    else
			    	return 0.0;
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

         
			inline void get_dense_matrix(DenseMatrix<T>& mat, const std::vector<Integer>& vec)const
			{
				Integer rows=vec.size();
				mat.resize(rows,rows);

				for(Integer i=0;i<rows;i++)
					for(Integer j=0;j<rows;j++)
                        {
                        	mat(i,j)=get_element_or_zero(vec[i],vec[j]);
                        	// std::cout<<"dense matrix i,j = "<<vec[i]<<", "<<vec[j]<<", mat="<<mat(i,j)<<std::endl;

                        }
				// for(Integer i=0;i<rows;i++)
					// {
					// for(Integer j=0;j<rows;j++)
                        
     //                    	std::cout<<mat(i,j)<<" ";
     //                 std::cout<<std::endl;

     //                 }


			}


			inline Integer max_cols_used()const {return max_cols_used_;}

			inline Integer max_rows()const {return max_rows_;}

			inline Integer max_cols()const {return max_cols_;}

			inline void set_max_cols_used(const Integer i) {max_cols_used_=i;}

			inline void set_max_rows(const Integer i) {max_rows_=i;}

			inline void set_max_cols(const Integer i) {max_cols_=i;}

			inline auto& cols_idx()const {return cols_idx_;}

			inline auto& rows_idx()const {return rows_idx_;}

			inline auto& cols_idx() {return cols_idx_;}

			inline auto& rows_idx() {return rows_idx_;}

			inline auto& non_zeros() {return non_zeros_;}




		};


	template<typename Mat,typename Vec,typename VecVec>
			void patch_gauss_seidel(Vec& x, const Mat& A, const Vec& b,
				                    const VecVec& entity2dofs,
				                    const Integer max_iter)
			{
                const int n=b.size();
                DenseMatrix<Real> local_mat;
                std::vector<Real> local_b;
                std::vector<Real> local_correction;
                std::vector<Real> local_rhs;
                Integer e2d_size;
                // todo fixme change initial conditions this works only for bc=0
				// if(x.size()!=n)
				// {
				// 	x.resize(n);
				// 	for(std::size_t i=0; i<n;i++)
				// 		x[i]=0;
				// }

				for(std::size_t it=0;it<max_iter;it++)
				{
					for(Integer i=0;i<entity2dofs.size();i++)
					{
						// std::cout<<"start i=="<<i<<std::endl;
						auto& e2d=entity2dofs[i];
						e2d_size=e2d.size();



						if(e2d_size)
						{
						// for(Integer j=0;j<e2d_size;j++)
						// 	std::cout<<e2d[j]<<" ";
						// std::cout<<std::endl;

						
						A.get_dense_matrix(local_mat,e2d);
                        // std::cout<<std::endl;
                        // std::cout<<"local_mat "<<std::endl;

         //                for(Integer j1=0;j1<local_mat.rows();j1++)
         //                {
						   // for(Integer j2=0;j2<local_mat.cols();j2++)
						   // 	std::cout<<local_mat(j1,j2)<<" ";
						   // std::cout<<std::endl;

         //                }

						subvector(local_b,x,e2d);
						A.multiply(local_rhs,x,e2d);


						// std::cout<<std::endl;

						// std::cout<<"e2d i=="<<i<<std::endl;
						// for(Integer j=0;j<e2d_size;j++)
						// 	std::cout<<e2d[j]<<" ";

						// std::cout<<std::endl;


						// std::cout<<"pre local_rhs i=="<<i<<std::endl;
						// for(Integer j=0;j<e2d_size;j++)
						// 	std::cout<<local_rhs[j]<<" ";

						// std::cout<<std::endl;

						// std::cout<<"b i=="<<i<<std::endl;
						// for(Integer j=0;j<e2d_size;j++)
						// 	std::cout<<b[e2d[j]]<<" ";

						// std::cout<<std::endl;

						for(Integer j=0;j<e2d_size;j++)
							local_rhs[j]=b[e2d[j]]-local_rhs[j];

						// std::cout<<std::endl;


						// std::cout<<"after local_rhs i=="<<i<<std::endl;
						// for(Integer j=0;j<e2d_size;j++)
						// 	std::cout<<local_rhs[j]<<" ";
                        


						// 	std::cout<<std::endl;


						local_correction.resize(e2d_size);
						local_mat.cholesky(local_correction,local_rhs);

						for(Integer j=0;j<e2d_size;j++)
							x[e2d[j]]=x[e2d[j]]+local_correction[j];

						// std::cout<<"local_correction i=="<<i<<std::endl;
						// for(Integer j=0;j<e2d_size;j++)
						// 	std::cout<<local_correction[j]<<" ";

						// std::cout<<"end i=="<<i<<std::endl;


						}





					}
					// std::cout<<"end iter=="<<it<<std::endl;


				}
				// std::cout<<"after iterations"<<std::endl;




			}



template<typename Elem,Integer ElementOrder, typename Operator,Integer FEFamily,Integer Order,Integer Continuity,Integer NComponents>
class Variable
{
public:  
// using Points=ElemPoints<Elem,Operator,FEFamily,Order>;
	using Points=ElemGeometricPoints<Elem,ElementOrder>;
	static constexpr Integer Npoints=Points::type::Dim;
	using RTnBaseFunctionSpace=BaseFunctionSpace<FEFamily,Order,Continuity, NComponents>;
	using BaseFunctionSpace=BaseElementFunctionSpace<Elem,FEFamily,Order,Continuity, NComponents>;
	using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;

    static constexpr Integer NDofs=FunctionSpaceDofsPerElem<FunctionSpace>::value;
    static constexpr Integer ManifoldDim=Elem::ManifoldDim;


	using SingleType   = typename SingleTypeShapeFunction<FunctionSpace,Operator>::SingleType;
	// using Coefficients=ShapeFunctionCoefficientsCollectionSingleType<Elem,FEFamily>;
	static constexpr Integer Rows=SingleType::Rows;
	static constexpr Integer Cols=SingleType::Cols;
	static constexpr Integer single_shape_function_components_number=Rows*Cols;
	static constexpr Integer solution_array_size=Npoints*Rows*Cols;
	using TotType = typename SingleTypeShapeFunction<FunctionSpace,Operator>::TotType;
	static constexpr Integer Ntot=FunctionSpaceDofsPerSubEntityElem<ElemFunctionSpace<Elem,BaseFunctionSpace>,Elem::ManifoldDim>::value;
	static constexpr Integer NDofs_single_space=Ntot/NComponents;
	static constexpr auto 
	reference_values{reference_shape_function_init<Elem,Operator,FEFamily,Order,SingleType,NDofs_single_space>(Points::points)};
	using Map=MapFromReference<Operator,Elem,FEFamily>;  
	static constexpr auto NShapes=reference_values.size();

    // using QRule=typename QuadratureRule<GaussianQuadrature>:: template rule<Elem,1>;
    using QRule=UserDefinedQuadratureRule<Elem,Points::type::Dim,ManifoldDim>;
    // using QRuleFace=Boundary2VolumetricQuadratureRule<QRuleFaceReference>;
    using Shape=ShapeFunction<Elem,RTnBaseFunctionSpace,IdentityOperator,QRule>;

    static constexpr auto MatPoints=vecmat2mat(Points::points);
    // using ShapeFaceValue=typename ShapeFace::type;
    // using QPpointsFace=typename ShapeFace::qp_points_type;


// qui devi passare i nodi dell'elemento
// l'ouput e' lungo elem:n_nods()*ncomponents

// EACH COMPONENT IS TREATED SEPARATELY, se vec is the array of dofs related to one component
        
        template<Integer Family>
		std::enable_if_t<Family==RaviartThomasFE, void> 
		init_aux( FiniteElem<Elem>& FE,Array<Real,Ntot> vec,const Integer component)
		{

		map_.init(FE);
		const auto& map=map_();
		// std::cout<<"component="<<component<<std::endl;
		// we compute u= sum_i u_i Map(phi_i)=Map(sum_i u_i phi_i) for each dof point on the element
		Integer cont=0;
        auto& mesh_ptr=FE.mesh_ptr();
        auto& signed_normal=mesh_ptr->signed_normal();
		SingleShapeFunctionCoefficientsCollection<Elem,FEFamily,Order>::apply(signed_normal.sign(FE.elem_id()),coeff_);
		shape_.init_map(map_);



		shape_.init(coeff_,MatPoints,FE);
		const auto& shape_eval=shape_.eval();
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
			
		// mat_tmp_[k](i,j)=vec[0+component]*reference_values[0][k](i,j);

		// mat_tmp_[k](i,j)=vec[0+component]*coeff_[0]*reference_values[0][k](i,j);

		mat_tmp_[k](i,j)=vec[0+component]*shape_eval[0][k](i,j);


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
				// mat_tmp_[k](i,j)+=vec[qp*NComponents+component]*reference_values[qp][k](i,j);

				// mat_tmp_[k](i,j)+=vec[qp*NComponents+component]*coeff_[qp]*reference_values[qp][k](i,j);


				mat_tmp_[k](i,j)+=vec[qp*NComponents+component]*shape_eval[qp][k](i,j);
			}
		// if(FE.elem_id()==18 ||FE.elem_id()==19 )

		// {
		// 	std::cout<<"FE.jac()"<<std::endl;
		// 	std::cout<<FE.jac()<<std::endl;
		// 	std::cout<<"FE.v0()"<<std::endl;
		// 	std::cout<<FE.v0()<<std::endl;
		// 	std::cout<<"coeff_"<<std::endl;
		// 	std::cout<<coeff_<<std::endl;
		// 	std::cout<<"reference_values"<<std::endl;
		// 	std::cout<<reference_values<<std::endl;
		// 	std::cout<<"vec"<<std::endl;
		// 	std::cout<<vec<<std::endl;
		// 	std::cout<<"pre mat_tmp_"<<std::endl;
		// 	std::cout<<mat_tmp_<<std::endl;
		// 	std::cout<<"map"<<std::endl;
		// 	std::cout<<map<<std::endl;
		// 	std::cout<<"shape_eval"<<std::endl;
		// 	std::cout<<shape_eval<<std::endl;
		// 	std::cout<<"MatPoints"<<std::endl;
		// 	std::cout<<MatPoints<<std::endl;
			
			
		// }

			// mat_tmp_[k]=map*mat_tmp_[k];




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

		// if(FE.elem_id()==18 ||FE.elem_id()==19 )

		// {
		// 	// std::cout<<"post mat_tmp_"<<std::endl;
		// 	// std::cout<<mat_tmp_<<std::endl;
		// 	// std::cout<<"output_"<<std::endl;
		// 	// std::cout<<output_<<std::endl;
		// }
		// std::cout<<"Npoints="<<Npoints<<std::endl;
		// std::cout<<"NComponents="<<NComponents<<std::endl;
		// std::cout<<"map="<<map<<std::endl;
		// std::cout<<"component="<<component<<std::endl;
		// std::cout<<"vec="<<vec<<std::endl;
		// std::cout<<"mat_tmp_="<<mat_tmp_<<std::endl;
		// std::cout<<" all output_="<<output_<<std::endl;
		// std::cout<<"reference_values="<<reference_values<<std::endl;
		}

        template<Integer Family>
		std::enable_if_t<Family==LagrangeFE, void> 
		init_aux( FiniteElem<Elem>& FE,Array<Real,Ntot> vec,const Integer component)
		{

		map_.init(FE);
		const auto& map=map_();
		Integer cont=0;
 
		for(Integer k=0;k<Npoints ;k++)
		{

		for(Integer i=0;i<Rows;i++)
		for(Integer j=0;j<Cols;j++)            
		{
		mat_tmp_[k](i,j)=vec[0+component]*reference_values[0][k](i,j);
		}


		for(Integer qp=1;qp<NShapes;qp++)
		for(Integer i=0;i<Rows;i++)
			for(Integer j=0;j<Cols;j++)            
			{
				mat_tmp_[k](i,j)+=vec[qp*NComponents+component]*reference_values[qp][k](i,j);
			}


			mat_tmp_[k]=map*mat_tmp_[k];


			for(Integer i=0;i<SingleType::Rows;i++)
				for(Integer j=0;j<SingleType::Cols;j++)
				{
					output_[cont]=mat_tmp_[k](i,j);
					cont++;
				}

			}
		}


        void init( FiniteElem<Elem>& FE,Array<Real,Ntot> vec,const Integer component)
        {
        	init_aux<FEFamily>(FE,vec,component);
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
		Array<Real,NDofs_single_space> coeff_;
		Shape shape_;
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
		// std::cout<<"print_solutions_aux_aux="<<std::endl;
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
				// std::cout<<"i=="<<i<<" n_elements=="<<n_elements<<std::endl;
                

				FE.init(i,level);
				dofsdofmap.template dofmap_get<N>(localdofmap,i,level);
                // std::cout<<"elem id=="<<i<<std::endl;
				// std::cout<<"localdofmap=="<<localdofmap<<std::endl;

				// std::cout<<"sol=="<<std::endl;
				// for(Integer k=0;k<sol.size();k++)
				// 	std::cout<<sol[k]<<std::endl;
				// std::cout<<"local_sol=="<<std::endl;
				// for(Integer k=0;k<local_sol.size();k++)
				// 	std::cout<<local_sol[k]<<std::endl;
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
		// std::cout<<"N="<<N <<std::endl;

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

    // std::cout<<"before print_solutions"<<std::endl;
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
	    		// std::cout<<"1 begin update N="<<N<<std::endl;
	    		node2elem_.init();
	    		node2elem_update_=true;
	    	}
	        if(!is_updated_[N]&&level_==-1)
	    	{
	    	// std::cout<<"2 real update N="<<N<<std::endl;
	    	is_updated_[N]=true;
	    	auto& ens=tuple_get<N>(tuple_entities_);
	    	ens.init_elem_entity(*mesh_ptr_,node2elem_.val(),level_);
	        }
	    	else if(!is_updated_[N]&&level_!=-1)
	    	{
	    	// std::cout<<"3 real update N="<<N<<std::endl;
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





	// template<typename MeshT_>
	// class MeshAndEntity
	// {
	//  public:
	//  	using MeshT=MeshT_;
	//  	using Elem=typename MeshT::Elem;
	//  	MeshAndEntity(const MeshT& mesh, const Integer level=-1):
	//     mesh_ptr_(std::make_shared<MeshT>(mesh)),
	//     level_(level),
	//     entities_(mesh,level)
	//  	{}


	//     inline auto mesh_ptr(){return mesh_ptr_;}

	//     inline auto level(){return level_;}

	//     inline auto bisection_ptr() {return entities_.bisection_ptr();}

	//     inline auto& entities_collection(){return entities_;}

	//     template<template<class...>class TemplateClass,typename...Args>
	//     inline void init(const TemplateClass<Args...>& W){return entities_.init(W);}

	//     template<typename Elem,typename FunctionSpace, typename...FunctionSpaces>
	//     inline void init(){return entities_.template init<Elem,FunctionSpace,FunctionSpaces...> ();}

	// 	inline void add_bisection(const Bisection<MeshT>& bisection) {entities_.add_bisection(bisection);}

		


	//  private:
	//  	std::shared_ptr<MeshT> mesh_ptr_;
	//  	Integer level_;
	//  	ElemEntityCollection<MeshT> entities_;
	 	
	// };



	// template<Integer N,typename MeshT>
	// auto& get_entity(MeshAndEntity<MeshT>& mesh_and_entity){return tuple_get<N>(mesh_and_entity.entities_collection().tuple_entities());}


	// template<typename Elem, typename...FunctionSpaces,typename MeshT>
	// void init(MeshAndEntity<MeshT>& mesh_and_entity){mesh_and_entity.template init<Elem,FunctionSpaces...>();}


















	// template<typename MeshT,Integer EntityDim>
	// class EntitySimplicialMap
	// {
	//  public:
	//  using IntegerVector=std::vector<Integer>;
	//  using Key=std::array<Integer,EntityDim+1>;
	//  using Value=std::shared_ptr<IntegerVector>;
	//  using Map=std::map<Key, Value>;
	//  using MapDof=std::map<Key, std::shared_ptr<Integer>>;
	//  static constexpr Integer Dim=MeshT::Dim;
	//  static constexpr Integer ManifoldDim=MeshT::ManifoldDim;
	//  static constexpr auto entity_combinations= ElemEntity<Simplex<Dim,ManifoldDim>,EntityDim>::entity_combinations();
	//  EntitySimplicialMap(const MeshT& mesh):
	//  mesh_ptr_(std::make_shared<MeshT>(mesh))
	//  {}

	//  void init()
	//  {
	//   Integer entity[EntityDim+1];
	//   Integer cont=0;
	//   for(std::size_t i=0;i<mesh_ptr_->n_elements();i++)
	// 	{
	// 		const auto& nodes=mesh_ptr_->elem(i).nodes;

	//           for(Integer iter_entity=0;iter_entity<entity_combinations;iter_entity++)
	//           {

	//             Combinations<ManifoldDim + 1, EntityDim+1>::generate(iter_entity,entity);

	//             for(int nn=0;nn<EntityDim+1;nn++)
	//               local_nodes_[nn]=nodes[entity[nn]];   

	//             std::sort(std::begin(local_nodes_), std::end(local_nodes_)); 


	//             auto& val=map_[local_nodes_];
	// 				if(!val) {
	// 					val = std::make_shared<IntegerVector>();
	// 				}
	//             auto& val2=entity_[local_nodes_];
	// 				if(!val2)
	// 				{

	// 					val2 = std::make_shared<Integer>(cont);
	// 					// entity_.insert(std::pair<Key,Integer>(local_nodes_,cont));
	// 					// entity_[local_nodes_]=cont;
	// 					// std::cout<<" local_nodes_="<<std::endl;
	// 					// for(int j=0;j<local_nodes_.size();j++)
	// 					// std::cout<<local_nodes_[j]<<" ";
	// 				 //    std::cout<<" CONT="<< entity_[local_nodes_]<<std::endl;

	// 					cont++;
	// 				}
	//             val->push_back(i);
	//             // .insert( std::pair<Key,Value>(local_nodes_,i)); 
	//         // for(std::size_t j=0;j<local_nodes_.size();j++)
	//         	// std::cout<<local_nodes_[j]<<" ";
	//         // std::cout<<std::endl;
	//         }


	// 	}

	//  }
	//  		void describe(std::ostream &os) const
	// 		{
	// 			os << "-----------------------------------\n";

	// 			for(const auto &m : map_) {
	// 				const auto& nodes=m.first;
	// 				for(const auto n : nodes) {
	// 					os << n << " ";
	// 				}

	// 				os << "-> ";

	// 				for(const auto n : *m.second) {
	// 					os << n << " ";
	// 				}

	// 				os << "\n";
	// 			}

	// 			os << "-----------------------------------\n";
	// 		}
	//  		void describe2(std::ostream &os) const
	// 		{
	// 			os << "----------------------------------- entity nums -------------------------\n";

	// 			for(const auto &m : entity_) {
	// 				const auto& nodes=m.first;
	// 				os << "NODES= ";
	// 				for(const auto n : nodes) {
	// 					os << n << " ";
	// 				}

	// 				os << "-> ";

	// 					os << *m.second << " ";
					

	// 				os << "\n";
	// 			}

	// 			os << "-----------------------------------\n";
	// 		}

	//  private:
	//  	std::shared_ptr<MeshT> mesh_ptr_;
	//  	Map map_;
	//  	std::array<Integer, EntityDim+1> local_nodes_;
	//     MapDof entity_;

	  
	// };




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


	 
	 Node2ElemMap(MeshT& mesh,BisectionT& bisection):
	 mesh_(mesh),
	 bisection_(bisection),
	 last_n_elements_(0)
	 {}
     

     Integer max_n_nodes()
     {
     	Integer n=0;
     	Integer n_size;
     	Integer cont;
     	for (auto it = map_.begin(); it != map_.end(); it++ )
			{   
				auto& second=it->second;
			    n_size=it->second->size();
                cont=0;
			    for(Integer i=0;i<n_size;i++)
			    	if(mesh_.is_active(second->operator[](i)))
			    		cont++;
				if(n<cont)
					n=cont;
			}
	     return n;

     }
	 // inline void add_bisection(const Bisection<MeshT>& bisection) {bisection_ptr_=std::make_shared<BisectionT>(bisection);}

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
	 		if(elem_belongs_to_level(mesh_,(*val)[i],level,bisection_.tracker()))
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
		BisectionT& bisection_;
		// std::shared_ptr<BisectionT> bisection_ptr_;
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
















	// template<typename MeshT_,Integer EntityDim>
	// class ConnectivitySimpliacialMap;

	// template<typename MeshT_,Integer EntityDim>
	// class ConnectivitySimpliacialMap
	// {

	// public:
	//  using MeshT=MeshT_;
	//  using BisectionT=Bisection<MeshT>;
	//  using Elem=typename MeshT::Elem;
	//  static constexpr Integer ManifoldDim=MeshT::ManifoldDim;
	//  static constexpr Integer Dim=MeshT::Dim;
	//  static constexpr auto entity_combinations= ElemEntity<Elem,EntityDim>::entity_combinations();
	//  static constexpr auto Npoints=EntityDim+1;
	//  using IntegerVector=std::vector<Integer>;
	//  using Key=std::array<Integer,Npoints>;
	//  using Value=std::shared_ptr<IntegerVector>;
	//  using Map=std::map<Key, Value>;
	//  using MapDof=std::map<Key,Integer>;

	 
	//  ConnectivitySimpliacialMap(MeshT& mesh,Node2ElemMap<MeshT>& node2elem,Bisection<MeshT>& bisection):
	//  mesh_(mesh),
	//  node2elem_(node2elem),
	//  bisection_(bisection),
	//  last_n_elements_(0),
	//  count_n_entity_(0),
	//  init_coarse_level_(false) 
	//  {
	//  	count_n_entity_vec_.resize(1,0);
	//  }




	//  void init_elem(const Elem elem,const Integer level )
	//  {

	//         // std::cout<<"___________LEVEL="<<level<<std::endl;
	//         // std::cout<<"count_n_entity_vec_.size()-1="<<count_n_entity_vec_.size()-1<<std::endl;
	//  		if(count_n_entity_vec_.size()-1<level)
	//  			count_n_entity_vec_.push_back(count_n_entity_vec_[count_n_entity_vec_.size()-1]);
	//  		// std::cout<<"___________count_n_entity_vec_[END]="<<std::endl;

	//  		// for(std::size_t i=0;i<count_n_entity_vec_.size();i++)
	//  			// std::cout<<count_n_entity_vec_[i]<<" ";
	//  		// std::cout<<std::endl;
	//  		Integer cont=0;
	//  		Integer cont1=0;
	//  		Integer cont2=0;

	//  	 	const auto& parent_id=elem.parent_id;
	// 	 	const auto& parent_elem=mesh_.elem(parent_id);
	// 	 	const auto& parent_nodes=parent_elem.nodes;
	// 	 	std::array<Integer,Npoints> parent_entity_nodes;
	// 	 	std::array<Integer,Npoints> child_entity_nodes;
	// 	 	std::array<Integer,Npoints> entity_used;
	//         bool found_parent_entity;
	 	 
	// 	 	 // find the parent element
	// 	 	 // loop on all its entities parent_entity of dimension EntityDim
	// 	 	 // given entity parent_entity:
	// 	 	 //     loop on all the child elements child_elem
	// 	 	 //          loop on all the entities child_entity
	// 	 	 //               if child_entity==parent_entity
	//          //                  do nothing, such entity already exists
	//          //      if no agreement has been found,  cont++
	//          // then we start creating new dofs for the entities on this level
	//          // entities which belong to both coarse and fine level are untouched
	//          // the new ones are created, counting them from n_coarse_entities-cont
	// 	 	 // we loop on all the entities of dimension EntityDim of the element
	// 	 	 const auto& child=parent_elem.children;

	//          for(std::size_t parent_entity=0;parent_entity<entity_combinations;parent_entity++)
	//          	{
	// 			 Combinations<ManifoldDim + 1, EntityDim+1>::generate(parent_entity,entity_);
	//              // std::cout<<std::endl;
	//          	 for(std::size_t i=0;i<Npoints;i++)
	//                 	{
	//                 		parent_entity_nodes[i]=parent_nodes[entity_[i]];
	//                 		// std::cout<<parent_entity_nodes[i]<<" ";
	//                 	}
	//              // std::cout<<std::endl;
	//              std::sort(parent_entity_nodes.begin(),parent_entity_nodes.end());
	            

	//             if(!parent_map_[parent_entity_nodes])
	//             {
	//             	// std::cout<<" entro in !parent_map_"<<std::endl;
	//             found_parent_entity=false;
	//             // loop on all the child elements child_elem          
	// 		 	for(std::size_t i=0;i<child.size();i++)
	// 		 	{
	// 		 	 const auto& child_elem=mesh_.elem(child[i]);
	// 		 	 const auto& child_id=child_elem.id;
	// 		 	 const auto& child_nodes=child_elem.nodes;
	//              // loop on all the entities child_entity
	// 	         for(std::size_t child_entity=0;child_entity<entity_combinations;child_entity++)
	// 	         	{
	// 				 Combinations<ManifoldDim + 1, EntityDim+1>::generate(child_entity,entity_);

	// 	         	 for(std::size_t i=0;i<Npoints;i++)
	// 	                	child_entity_nodes[i]=child_nodes[entity_[i]];
	// 	             std::sort(child_entity_nodes.begin(),child_entity_nodes.end());
	//                  // parent_map_[parent_entity_nodes]=true;
	// 	             found_parent_entity=std::equal(child_entity_nodes.begin(),child_entity_nodes.end(),parent_entity_nodes.begin());
				     
	// 			     if(found_parent_entity)			     	 
	// 			     	 goto label;
				     	
	// 			     }
	// 			  }
	// 			  label:
	// 			  if(found_parent_entity)
	// 			  {
	// 			  	// parent_map_[parent_entity_nodes]=true;
	// 			  	// std::cout<<"found_parent_entity="<<found_parent_entity<<std::endl;
	// 			  	// entity_used[cont1]=map_dof_[parent_entity_nodes];
	// 			  	// cont1++;
	// 			  }
	// 			  else
	// 			  {
	// 			     	// std::cout<<"entity_used["<<cont1<<"]="<<map_dof_[parent_entity_nodes]<<std::endl;
	//                     parent_map_[parent_entity_nodes]=true;
	// 			  	    entity_used[cont]=map_dof_[parent_entity_nodes];
	// 			  		cont++;
	// 			  }
	// 			  {}
	// 			  }
	// 	         }


	//             count_n_entity_vec_[level]=count_n_entity_vec_[level];//-cont;
	//             // loop on all the child elements child_elem          
	// 		 	for(std::size_t i=0;i<child.size();i++)
	// 		 	{
	// 		 	 const auto& child_elem=mesh_.elem(child[i]);
	// 		 	 const auto& child_id=child_elem.id;
	// 		 	 const auto& child_nodes=child_elem.nodes;			 		
	// 		 	 if(elem2dofs_[child_id].size()==0)
	// 		 	 {

	//              // loop on all the entities child_entity
	// 	         for(std::size_t child_entity=0;child_entity<entity_combinations;child_entity++)
	// 	         	{
	// 				 Combinations<ManifoldDim + 1, EntityDim+1>::generate(child_entity,entity_);
	// 	         	 for(std::size_t i=0;i<Npoints;i++)
	// 	                	child_entity_nodes[i]=child_nodes[entity_[i]];
	// 	             std::sort(child_entity_nodes.begin(),child_entity_nodes.end());
	                 
	//                  auto& child_map=map_[child_entity_nodes];

	//                  // if the entity has not been visited yet
	//                  if(!child_map)
	//                  {
	//                  	// create the vector
	//                     child_map=std::make_shared<IntegerVector>();
	//                     auto& child_dof=map_dof_[child_entity_nodes];
	//                     if(cont2<cont)
	//                     {
	//                      child_dof=entity_used[cont2];
	//                      // std::cout<<" cont2<cont "<<child_dof<<std::endl;
	//                      cont2++;	
	//                     }
	//                     else
	//                     {child_dof=count_n_entity_vec_[level];
	//                      // std::cout<<" cont2>=cont "<<child_dof<<std::endl;

	// 					 count_n_entity_vec_[level]++;}
	//                  }
	//                  // std::cout<<"count_n_entity_vec_[level]="<<count_n_entity_vec_[level]<<std::endl;
	//                  child_map->push_back(child_id);
	//                  std::cout<<" child_id = "<<child_id<<" with "<<map_dof_[child_entity_nodes]<<std::endl;
	//                  elem2dofs_[child_id].push_back(map_dof_[child_entity_nodes]);
	// 			  }
	// 		 	 }
	// 			}
	// 	 	 }


	//  void init_coarse()
	//  {  
	//  	node2elem_.init();
	//  	elem2dofs_.resize(mesh_.n_elements());
	//  	for(std::size_t i=last_n_elements_;i<mesh_.n_elements();i++)
	//  	{
	//  		const auto& elem=mesh_.elem(i);
	//  		if(elem.parent_id==INVALID_INDEX)
	//  		{
	//  		last_n_elements_++;
	// 		const auto& nodes=mesh_.elem(i).nodes;

	//           for(Integer iter_entity=0;iter_entity<entity_combinations;iter_entity++)
	//           {
	//            // std::cout<<" ini coarse = "<<i<<std::endl;
	//             Combinations<ManifoldDim + 1, Npoints>::generate(iter_entity,entity_);

	//             for(int nn=0;nn<Npoints;nn++)
	//               local_nodes_[nn]=nodes[entity_[nn]];   

	//             std::sort(std::begin(local_nodes_), std::end(local_nodes_)); 


	//             auto& val=map_[local_nodes_];
	// 				if(!val) {
	// 					val = std::make_shared<IntegerVector>();
	// 					// parent_map_[local_nodes_]=true;
	// 				}
	//             auto& val2=map_dof_[local_nodes_];
	// 				if(!val2)
	// 				{

	// 					val2 = count_n_entity_;
	// 					count_n_entity_++;
	// 				}
	//             val->push_back(i);
	//             // std::cout<<"qui elem2dofs_.size="<<elem2dofs_.size()<<std::endl;
	//              // std::cout<<"elem.parent_id="<<elem.id<<std::endl;
 //                elem2dofs_[elem.id].push_back(val2);
 //                // std::cout<<" after qui elem2dofs_.size="<<elem2dofs_.size()<<std::endl;
	//           }
	//  		}
	//       }
	//      // std::cout<<"end1 ini coarse = "<<std::endl;
	//      count_n_entity_vec_[0]=count_n_entity_;
	//      // std::cout<<" end2 coarse = "<<std::endl;
	//      // count_levels_=1;
	//      // std::cout<<" end3 coarse = "<<std::endl;
	// 	}             

	//  void init()
	//  {
	//  	node2elem_.init();
	//  	elem2dofs_.reserve(mesh_.n_elements());

 //        // we initialize the coarse level only the first time
	//  	if(!init_coarse_level_)
	//  	{
	//  		init_coarse();
	//  		init_coarse_level_=true;
	//  	}
        
 //        // std::cout<<"elem2dofs_ with size="<<elem2dofs_.size()<<std::endl;
 //        // for(std::size_t i=0;i<elem2dofs_.size();i++)
 //        // {
 //        // 	std::cout<<std::endl;
 //        // 	for(std::size_t j=0;j<elem2dofs_[i].size();j++)
 //        // 	std::cout<<elem2dofs_[i][j]<<" ";
 //        // std::cout<<std::endl;
 //        // }



	//     // we assume the finest elements are added after the coarsest ones  
	//  	const auto& tracker=bisection_.tracker();
	//     Integer level;
	//     // std::cout<<"last_n_elements_="<<last_n_elements_<<std::endl;
	//     // std::cout<<"mesh_.n_elements()="<<mesh_.n_elements()<<std::endl;    
	//  	for(std::size_t i=last_n_elements_;i<mesh_.n_elements();i++)
	//  	   {
	 	   	
	//  	   	level=tracker.get_iterate(i);
	//  	   	// std::cout<<"init elem="<<i <<" of level "<<level<<std::endl;
	//  		init_elem(mesh_.elem(i),tracker.get_iterate(i));
	//  		// std::cout<<count_n_entity_vec_[level]<<std::endl;
	//  		last_n_elements_++;	
	//         }
	//     last_n_elements_=mesh_.n_elements();

	 	

	//   } 


	//  		void describe(std::ostream &os) const
	// 		{
	// 			os << "-----------------------------------\n";

	// 			for(const auto &m : map_) {
	// 				const auto& nodes=m.first;
	// 				for(const auto n : nodes) {
	// 					os << n << " ";
	// 				}

	// 				os << "-> ";

	// 				for(const auto n : *m.second) {
	// 					os << n << " ";
	// 				}

	// 				os << "\n";
	// 			}

	// 			os << "-----------------------------------\n";
	// 		}

	//  		void describe2(std::ostream &os) const
	// 		{
	// 			os << "-----------------------------------\n";

	// 			for(const auto &m : map_dof_) {
	// 				const auto& nodes=m.first;
	// 				for(const auto n : nodes) {
	// 					os << n << " ";
	// 				}

	// 				os << "-> ";

	// 					os << m.second << " ";
					

	// 				os << "\n";
	// 			}

	// 			os << "-----------------------------------\n";
	// 		}		



	//  		void describe3(std::ostream &os) const
	// 		{
	// 			os << "-----------------------------------\n";

	// 			for(const auto &m : parent_map_) {
	// 				const auto& nodes=m.first;
	// 				for(const auto n : nodes) {
	// 					os << n << " ";
	// 				}

	// 				os << "-> ";

	// 					os << m.second << " ";
					

	// 				os << "\n";
	// 			}

	// 			os << "-----------------------------------\n";
	// 		}	



	// 		auto& elem2dofs(const Integer i)const{return elem2dofs_[i];}


	// private:
	// 	MeshT& mesh_;
	// 	Node2ElemMap<MeshT>& node2elem_;
	// 	Bisection<MeshT>& bisection_;
	// 	std::shared_ptr<BisectionT> bisection_ptr_;
	//     Map map_;
	//     MapDof map_dof_;
	//     std::map<Key, bool> parent_map_;
	//     Key local_nodes_;
	//     Integer entity_[Npoints];
	//     Simplex<Dim,EntityDim> simplex_;
	//     Simplex<Dim,EntityDim> simplex_parent_;
	//     Integer last_n_elements_;
	//     Integer count_n_entity_;
	//     IntegerVector count_n_entity_vec_;
	//     // Integer count_levels_;
	//     bool init_coarse_level_;
	//     std::vector<std::vector<Integer>> elem2dofs_;
	// };



















	// template<typename Elem, Integer EntityDimMax, Integer EntityDim>
	//  class ConnectivitySimpliacialMapOfMeshTupleTypeHelper;

	// template<typename MeshT, Integer EntityDimMax, Integer EntityDim>
	// class ConnectivitySimpliacialMapOfMeshTupleTypeHelper
	// {
	// public:
	//      using rest = typename ConnectivitySimpliacialMapOfMeshTupleTypeHelper<MeshT,EntityDimMax,EntityDim+1>::type;
	//      using ens  = ConnectivitySimpliacialMap<MeshT,EntityDim>;
	//      using type = decltype( std::tuple_cat( std::declval< std::tuple<ens> >(),std::declval< rest >() ) );
	// };


	// template<typename MeshT, Integer EntityDimMax>
	// class ConnectivitySimpliacialMapOfMeshTupleTypeHelper<MeshT, EntityDimMax,EntityDimMax>
	// {
	// public:
	//      using ens  = ConnectivitySimpliacialMap<MeshT,EntityDimMax>;
	//      using type = typename std::tuple<ens>;
	// };

	// template<typename MeshT>
	// using ConnectivitySimpliacialMapOfMesh= typename ConnectivitySimpliacialMapOfMeshTupleTypeHelper<MeshT,MeshT::Elem::ManifoldDim,0>::type;



	// 		template<typename MeshT_>
	// class ConnectivitySimpliacialMapCollection
	// {
	// public:
	// 	using MeshT=MeshT_;
	// 	using BisectionT=Bisection<MeshT>;
	// 	using Elem=typename MeshT::Elem;
	// 	using EntitiesTuple=ConnectivitySimpliacialMapOfMesh<MeshT>;
	// 	static constexpr Integer ManifoldDim=Elem::ManifoldDim;


	// 		template<Integer N=0>
	// 	struct ConstructorTupleHelper;

	// 		template<>
	// 	struct ConstructorTupleHelper<0>
	// 	{
	// 		using ens = ConnectivitySimpliacialMap<MeshT,0>;
	// 		using type = typename std::tuple<ens>;
	// 	};

	// 		template <Integer N>
	// 	struct ConstructorTupleHelper
	// 	{
	// 		using rest = typename ConstructorTupleHelper<N-1>::type; 
	// 		using ens = ConnectivitySimpliacialMap<MeshT,N>;
	// 		using tuple_ens=std::tuple<ens>;
	// 		using type = decltype( std::tuple_cat(std::declval< rest >() , std::declval< tuple_ens >()) );
	// 	};

	// 		template<Integer N>
	// 	using ConstructorTuple=typename ConstructorTupleHelper<N>::type;




	// 		template< Integer N>
	// 	std::enable_if_t<0==N, ConstructorTuple<N> >
	// 	construct_tuple(MeshT& mesh,Node2ElemMap<MeshT>& node2elem,Bisection<MeshT>& bisection)
	// 	{   
	// 		using type=ConnectivitySimpliacialMap<MeshT,N>;
	// 		return std::tuple<type>
	// 		(type(mesh,node2elem,bisection));
	// 	}


	// 		template<Integer N>
	// 	std::enable_if_t< 0<N, ConstructorTuple<N> >
	// 	construct_tuple(MeshT& mesh,Node2ElemMap<MeshT>& node2elem,Bisection<MeshT>& bisection)
	// 	{
	// 		using type=ConnectivitySimpliacialMap<MeshT,N>;
	// 		return std::tuple_cat(construct_tuple<N-1>(mesh,node2elem,bisection), 
	// 			std::tuple<type>(type(mesh,node2elem,bisection)));
	// 	}








	// 	ConnectivitySimpliacialMapCollection(MeshT& mesh,Node2ElemMap<MeshT>& node2elem,Bisection<MeshT>& bisection):
	// 	mesh_(mesh),
	// 	node2elem_(node2elem),
	// 	bisection_(bisection),
	// 	tuple_entities_(construct_tuple<MeshT::ManifoldDim>(mesh,node2elem,bisection))
	// 	{}



	// 		    template<typename FunctionSpace,Integer Nmax,Integer N>
	// 	std::enable_if_t<(N>Nmax),void>
	// 	init_aux_aux()
	// 	{}

	// 		    template<typename FunctionSpace,Integer Nmax,Integer N>
	// 	std::enable_if_t<(N<=Nmax),void>
	// 	init_aux_aux()
	// 	{   
	// 		auto& ens=tuple_get<FunctionSpace::entity[N]>(tuple_entities_);
	// 		ens.init();
	// 		init_aux_aux<FunctionSpace,Nmax,N+1>();
	// 	}


	// 		    template<typename BaseFunctionSpace, typename...BaseFunctionSpaces>
	// 	std::enable_if_t<(sizeof...(BaseFunctionSpaces)==0),void>
	// 	init_aux()
	// 	{
	// 		using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;
	// 		init_aux_aux<FunctionSpace,FunctionSpace::entity.size()-1,0>();
	// 	}

	// 		    template<typename BaseFunctionSpace, typename...BaseFunctionSpaces>
	// 	std::enable_if_t<(sizeof...(BaseFunctionSpaces)>0),void>
	// 	init_aux()
	// 	{   
	// 		using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;
	// 		init_aux_aux<FunctionSpace,FunctionSpace::entity.size()-1,0>();
	// 		init_aux<BaseFunctionSpaces...>();
	// 	}


	// 		    template<typename BaseFunctionSpace, typename...BaseFunctionSpaces>
	// 	void init()
	// 	{   
	// 		init_aux<BaseFunctionSpace,BaseFunctionSpaces...>();
	// 	}


 //        template<Integer Nmax,Integer N>
	//  	std::enable_if_t<(N>Nmax),void> describe2_aux(std::ostream &os) const
	//  	{}

 //        template<Integer Nmax,Integer N>
	//  	std::enable_if_t<N<=Nmax,void> describe2_aux(std::ostream &os) const
	// 	{
	// 		os << "-----------------ENTITY N=="<<N<<"\n";
	// 		tuple_get<N>(tuple_entities_).describe2(os);
	// 		describe2_aux<Nmax,N+1>(os);
	// 	}

	//  	void describe2(std::ostream &os) const
	// 	{
 //         describe2_aux<ManifoldDim,0>(os);
	// 	}	


	// 	auto& tuple_entities(){return tuple_entities_;}

	// 	auto& node2elem(){return node2elem_;}

	// 	auto& mesh(){return mesh_;}

	// 	auto& bisection(){return bisection_;}

 //        const auto n_levels(){return bisection_.tracker().current_iterate();}
	// 	// auto& mesh_ptr(){return std::make_shared<MeshT>(mesh_);}

	// 	// inline auto bisection_ptr(){return bisection_ptr_;}
	// private:
	// 	MeshT& mesh_;
	// 	Node2ElemMap<MeshT>& node2elem_;
	// 	Bisection<MeshT>& bisection_;
	// 	// std::shared_ptr<BisectionT> bisection_ptr_;


	// 	bool node2elem_update_;
	// 	EntitiesTuple tuple_entities_;
	// 	// Array<bool, ManifoldDim+1> is_updated_;
	// 	Integer level_;
	// };






 //    template<typename...BaseFunctionSpaces,typename MeshT>
 //    void init(ConnectivitySimpliacialMapCollection<MeshT>& c) 
 //    {c.template init<BaseFunctionSpaces...>();}




template<typename ElemFunctionSpace>
class LagrangeDofs;


template<Integer Dim,Integer ManifoldDim,Integer Order, Integer Continuity, Integer NComponents>
class LagrangeDofs<ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,Order,Continuity,NComponents>>
{
public:
	using FunctionSpace=ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,Order,Continuity,NComponents>;
    static constexpr Integer NDofs=FunctionSpaceDofsPerElem<FunctionSpace>::value;
	static constexpr auto value()
	{
     const auto& ref=Simplex<Dim,ManifoldDim>::reference_init();
     Array<Matrix<Real, 1, Dim>,NDofs> dofs;
     const Integer Mref=ref.size();
     const Integer Nref=ref[0].size();
     Integer cont=Mref;
     for(Integer i=0;i<Mref;i++)
     	for(Integer j=0;j<Nref;j++)
        {
        	dofs[i](0,j)=ref[i][j];

        }

     if(Order==2)
     {
     for(Integer i=0;i<Mref;i++)
     	for(Integer j=i+1;j<Mref;j++)
     	{
	     	for(Integer k=0;k<Nref;k++)
	        {
	        	dofs[cont](0,k)=0.5*(ref[i][k]+ref[j][k]);
	        }

	        cont++;

     	}	
     }

      return dofs;



	}
};

template<typename ElementFunctionSpaceSpace,Integer QR=GaussianQuadrature>
class Dof;

template<typename ElementFunctionSpaceSpace,Integer QR=GaussianQuadrature>
class DofAux;




template <Integer Dim,Integer ManifoldDim, Integer Order,Integer Continuity, Integer NComponents_,Integer QR>
class DofAux<ElementFunctionSpace<Simplex<Dim,ManifoldDim>,LagrangeFE,Order,Continuity,NComponents_>,QR>
{
 public:

    using Elem=Simplex<Dim,ManifoldDim>;
 	using QuadratureElem=typename VolumeOrSurfaceElem<Elem,false>::type;
 	static constexpr Integer FEFamily=LagrangeFE;
 	static constexpr Integer NComponents=1; 
 	using BaseFunctionSpace=BaseFunctionSpace<FEFamily,Order,Continuity,NComponents>;

    static constexpr Integer QRuleOrder=QuadratureOrder<IdentityOperator, BaseFunctionSpace>::value;


    using QRule=typename QuadratureRule<QR>:: template rule<Elem,1>;
    // using Shape=ShapeFunction<Elem,BaseFunctionSpace,IdentityOperator,QRule>;
    using Shape=ShapeFunction<Elem,BaseFunctionSpace,IdentityOperator,QRule>;
    using ShapeValue=typename Shape::type;






 	// HA SENSO PER NCOMPONENTS == 1 RIFLETTI PER PIU COMPONENTI


 	using ElementFunctionSpace=ElementFunctionSpace<Elem,FEFamily,Order,Continuity,NComponents>;
 	static constexpr auto entity=ElementFunctionSpace::entity;
    using Map=MapFromReference<IdentityOperator,Elem,FEFamily> ;
    using MapReference=MapFromReference<IdentityOperator,QuadratureElem,FEFamily> ;
    static constexpr auto Dofs=LagrangeDofs<ElementFunctionSpace>::value();
    static constexpr auto NDofs=LagrangeDofs<ElementFunctionSpace>::NDofs;
    


    template<Integer M,Integer N=0,typename Mat>
    inline std::enable_if_t<N==NDofs,void>
    loop_coarse_dof(Mat& mat)
    {}


    template<Integer M,Integer N=0,typename Mat>
    inline std::enable_if_t<N<NDofs,void>
    loop_coarse_dof(Mat& mat)
    {
     // HA SENSO PER NCOMPONENTS == 1 RIFLETTI PER PIU COMPONENTI
     // TODO FIXME
     // auto& value=func_values_[N];

     // DofAux<M,N,ElementFunctionSpace,QR>::compute(local_mat_,func_values_);
     // std::cout<<"loop_coarse_dof =="<<N<<std::endl;
     // std::cout<<value<<std::endl;
     mat(M,N)=func_values_[N][0](0,0);

     // std::cout<<"(M,N) ==("<<M<<", "<<N<<")=>"<<func_values_[N][0](0,0)<<std::endl;
     // shape_.init();
     // loop on coarse shape functions
     loop_coarse_dof<M,N+1>(mat);
        
    }





    
    void compute_point_row(FiniteElem<Elem>& C_FE, FiniteElem<Elem>& F_FE)
    {
    const auto& J_inv=C_FE.inv_jac();
    const auto& v0=C_FE.v0();
    const auto& rows=J_inv.rows();
    const auto& cols=J_inv.cols();
    // std::cout<<J_inv<<std::endl;
    // std::cout<<v0<<std::endl;
    // std::cout<<rows<<std::endl;
    // std::cout<<cols<<std::endl;
    // std::cout<<"transformed_point_col_"<<std::endl;
    // std::cout<<transformed_point_col_<<std::endl;

    for(Integer i=0;i<rows;i++)
    	{
    		// std::cout<<"i="<<i<<std::endl;
    		transformed_point_row_(0,i)=J_inv(i,0)*(transformed_point_col_[0]-v0[0]);
    	}

    for(Integer i=0;i<rows;i++)
    	{
    		// std::cout<<"i="<<i<<std::endl;
	    	for(Integer j=1;j<cols;j++)
	    	{
	    		// std::cout<<"j="<<j<<std::endl;
	    		transformed_point_row_(0,i)+=J_inv(i,j)*(transformed_point_col_[j]-v0[j]);
	    	}
       	}

    // std::cout<<"transformed_point_row_"<<std::endl;
    // std::cout<<transformed_point_row_<<std::endl;
    //  std::cout<<"compute_point_row end"<<std::endl;

     }



    template<Integer N=0,typename Mat>
    inline std::enable_if_t<N==NDofs,void>
    loop_fine_dof(Mat& mat, FiniteElem<Elem>& C_FE, FiniteElem<Elem>& F_FE)
    {}

    template<Integer N=0,typename Mat>
    inline std::enable_if_t<N<NDofs,void>
    loop_fine_dof(Mat& mat, FiniteElem<Elem>& C_FE, FiniteElem<Elem>& F_FE)
    {
     // std::cout<<"loop_fine_dof =="<<N<<"    "<<Dofs[N]<<std::endl;
     // F_FE.init(F_FE.elem_id());
     F_map_.init(F_FE);
     F_FE.transform_point(transformed_point_col_,Dofs[N]);
     // std::cout<<"transformed_point_col_"<<std::endl;
     // std::cout<<transformed_point_col_<<std::endl;
     // transformed_point_row_=C_FE_.inv_jac()*transformed_point_col_;
     compute_point_row(C_FE,F_FE);
     // std::cout<<"shape init first"<<std::endl;

     shape_.init_map(C_map_);
     shape_.init(transformed_point_row_,C_FE);
     // std::cout<<"shape init end"<<std::endl;




     // std::cout<<"Dofs[N] =="<<std::endl;
     // std::cout<<Dofs[N]<<std::endl;
     // std::cout<<"transformed_point_col_ =="<<std::endl;
     // std::cout<<transformed_point_col_ <<std::endl;
     // std::cout<<"transformed_point_row_ =="<<std::endl;
     // std::cout<<transformed_point_row_ <<std::endl;
     func_values_=shape_.eval();
     // std::cout<<"func_values_ =="<<std::endl;
     // std::cout<<func_values_<<std::endl;

     loop_coarse_dof<N>(mat);


     loop_fine_dof<N+1>(mat, C_FE,F_FE);
        
    }

    // template<typename Mat>
    // inline void compute(Mat&mat, FiniteElem<Elem>& C_FE, FiniteElem<Elem>& F_FE)
    // {

    //  // C_FE_.init(C_FE.elem_id());
    //  C_FE.init_inv_jac();
    //  C_map_.init(C_FE);
    //  shape_.init_map(C_reference_map_);



    //  // loop on coarse shape functions
    //  // std::cout<<"compute"<<std::endl;
    //  // std::cout<<"map_"<<C_map_()<<std::endl;
    //  // std::cout<<"Dofs reference =="<<std::endl;
    //  // std::cout<<Dofs<<std::endl;
    //  loop_fine_dof(mat,C_FE,F_FE);
    //  std::cout<<"local_mat_"<<std::endl;
    //  std::cout<<local_mat_<<std::endl;
    // }








    template<typename Mat>
    void compute(Mat& mat, FiniteElem<Elem>& C_FE, FiniteElem<Elem>& F_FE)
    {
    	Integer cont=0;

    	C_FE.init_inv_jac();
    	C_map_.init(C_FE);
    	// compute_aux<Cont>(mat,cont,C_FE,F_FE);
    	// compute_face_aux<0>(mat,cont,C_FE,F_FE);
    	// compute_momentum_aux<0>(mat,cont,C_FE,F_FE);

    	loop_fine_dof(mat,C_FE,F_FE);


        // std::cout<<"mat"<<std::endl;
        // std::cout<<mat<<std::endl;
    	// for(Integer i=0;i<ManifoldDim+1;i++)
    	// 	{F_FE.init_boundary(i,true);
    	//      compute(mat,cont,C_FE,F_FE);}


    	// compute<entity[1]>(mat,C_FE,F_FE);
    }  

  // template<typename Mat,typename Values>
  // static void compute(Mat& local_mat, Values& func_values )
  // {
  // 	local_mat(M_,N_)=func_values[N_][0](0,0);
  // }

  private:
 	Shape shape_;
 	Map C_map_;
 	Map F_map_;
 	// MapReference C_reference_map_;
    ShapeValue func_values_;
    Matrix<Real,1,Dim> transformed_point_row_;
    Matrix<Real,Dim,1> transformed_point_col_;
  //   Matrix<Real,NDofs,NDofs> local_mat_;
};


template <Integer Dim,Integer ManifoldDim,Integer Continuity, Integer NComponents_,Integer QR>
class DofAux<ElementFunctionSpace<Simplex<Dim,ManifoldDim>,RaviartThomasFE,0,Continuity,NComponents_>,QR>
{
 public:

    using Elem=Simplex<Dim,ManifoldDim>;
    static constexpr Integer Order=0;
 	using Face=typename VolumeOrSurfaceElem<Elem,false>::type;
 	static constexpr Integer FEFamily=RaviartThomasFE;
 	static constexpr Integer NComponents=1; 
 	using RTnBaseFunctionSpace=BaseFunctionSpace<FEFamily,Order,Continuity,NComponents>;
 	using ElementFunctionSpace=ElementFunctionSpace<Elem,FEFamily,Order,Continuity,NComponents>;
 	static constexpr auto entity=ElementFunctionSpace::entity;
 	static constexpr auto dofs_per_entity=ElementFunctionSpace::dofs_per_entity;

    static constexpr Integer NDofs=FunctionSpaceDofsPerElem<ElementFunctionSpace>::value;

    static constexpr Integer NDofsFace= NComponents * dofs_per_entity[0] * ElemEntityCombinations<Elem,entity[0]>::value;
    static constexpr Integer NDofsMomentum= NComponents * dofs_per_entity[1] * ElemEntityCombinations<Elem,entity[1]>::value;

    // we have to integrate so + 1
    static constexpr Integer QRuleOrder=1;//QuadratureOrder<IdentityOperator, RTnBaseFunctionSpace>::value+1;
    

    using Map=MapFromReference<IdentityOperator,Elem,FEFamily> ;
    using TraceMap=MapFromReference<TraceOperator,Face,FEFamily> ;

    using F_Lagrangian_Map=MapFromReference<IdentityOperator,Elem,LagrangeFE> ;


    using QRuleFaceReference=typename QuadratureRule<QR>:: template rule<Face,QRuleOrder>;
    using QRuleFace=Boundary2VolumetricQuadratureRule<QRuleFaceReference>;
    using ShapeFace=ShapeFunction<Elem,RTnBaseFunctionSpace,IdentityOperator,QRuleFace>;
    using ShapeFaceValue=typename ShapeFace::type;
    using QPpointsFace=typename ShapeFace::qp_points_type;


    using ShapeTraceFace=ShapeFunction<Elem,RTnBaseFunctionSpace,TraceOperator,QRuleFaceReference>;
    using trace_type=typename ShapeTraceFace::trace_type;

    // using QRuleMomentum=typename QuadratureRule<QR>:: template rule<Elem,QRuleOrder>;
    // using ShapeMomentum=ShapeFunction<Elem,RTnBaseFunctionSpace,IdentityOperator,QRuleMomentum>;
    // using ShapeMomentumValue=typename ShapeMomentum::type;
    // using QPpointsMomentum=typename ShapeMomentum::qp_points_type;

    using Coefficients=ShapeFunctionCoefficientsCollectionSingleType<Elem,FEFamily>;
    // using ShapeLagrangianMomentum=ShapeFunction<Elem,PnBaseFunctionSpace,IdentityOperator,QRuleMomentum>;

  // template<typename Mat,typename Values>
  // static void compute(Mat& local_mat, Values& func_values )
  // {
  // 	local_mat(M,N)=func_values[N][0](0,0);
  // }
    template<typename Output, typename Weight, typename ShapeValue,typename Normal>
    void compute_face_dofs(Output& out, Integer& fine_dof, const Weight& weights,const ShapeValue& shape_value, const Normal& normal, FiniteElem<Elem>& FE)
    {
     const auto& NQPoints=ShapeValue::NQPoints;
     const auto& Ndofs=Output::Rows;

     // std::cout<<"compute_face_dofs"<<std::endl;
     // std::cout<<"shape_value"<<std::endl;
     // std::cout<<shape_value<<std::endl; 

     // std::cout<<"normal"<<std::endl;
     // std::cout<<normal<<std::endl;
     // std::cout<<"weights"<<std::endl;
     // std::cout<<weights<<std::endl;
     // std::cout<<"FE.get_det_side()"<<std::endl;
     // std::cout<<FE.get_det_side()<<std::endl;
     for(Integer coarse_dof=0;coarse_dof<NDofs;coarse_dof++)
     {
        out(fine_dof,coarse_dof)=0;

     	for(Integer qp=0;qp<NQPoints;qp++)
     	{
     		for(Integer i=0;i<ManifoldDim;i++)
     		{
     			out(fine_dof,coarse_dof)+= weights[qp]*shape_value[coarse_dof][qp][i]*normal[i];
     		}
     	}
       out(fine_dof,coarse_dof)*=FE.get_det_side();

     }
     // std::cout<<"out"<<std::endl;
     // std::cout<<out<<std::endl;    
    }

    // face dofs
    template<Integer EntityDim, typename Mat>
    std::enable_if_t< (EntityDim==entity[0]), void>
    compute(Mat& mat, Integer& cont, FiniteElem<Elem>& C_FE, FiniteElem<Elem>& F_FE)
    {
     // F_map_.init(F_FE_);
    	// QPpointsFace frfrfr(5);
    	// typename QRuleFaceReference::qp_points_type feef(5);
    	// std::cout<<"compute "<<std::endl;
    	auto mesh_ptr=F_FE.mesh_ptr();
    	// std::cout<<"qui1 "<<std::endl;
    	const auto& signed_normal=mesh_ptr->signed_normal();
    	// std::cout<<"qui2 "<<std::endl;
        C_map_.init(C_FE);
        // std::cout<<"qui3 "<<std::endl;
        shape_face_.init_map(C_map_);
        // std::cout<<"qui4 "<<std::endl;
        // coeffs_.init(*mesh_ptr);

        SingleShapeFunctionCoefficientsCollection<Elem,FEFamily,Order>::apply(signed_normal.sign(C_FE.elem_id()),coeff_);
        // std::cout<<"qui5 "<<std::endl;
    	F_FE.side_transform_point(qp_points_face_,QRuleFaceReference::qp_points);
    	// std::cout<<"qui6 "<<std::endl;
     // 	std::cout<<"QRuleFaceReference::qp_points"<<std::endl;
    	// std::cout<<QRuleFaceReference::qp_points<<std::endl;
    	// std::cout<<"qp_points_face_"<<std::endl;
    	// std::cout<<qp_points_face_<<std::endl;
     // transformed_point_row_=C_FE_.inv_jac()*transformed_point_col_;
        C_FE.points_to_reference_points(reference_qp_points_face_,qp_points_face_);
        // std::cout<<"qui7 "<<std::endl;
    	// std::cout<<"reference_qp_points_face_"<<std::endl;
    	// std::cout<<reference_qp_points_face_<<std::endl;
        shape_face_.init(coeff_,reference_qp_points_face_,C_FE);
        // std::cout<<"shape_face_.eval()"<<std::endl;
        // std::cout<<shape_face_.eval()<<std::endl;
        // // Vector<Real,NDofsFace> out; 

        const auto& normal=signed_normal.normals()[F_FE.elem_id()][F_FE.side_id()];
        

        // std::cout<<"shape_trace_face_.eval()"<<std::endl;
        // std::cout<<shape_trace_face_.eval()<<std::endl;
        
         // std::cout<< "F_FE elem id = "<< F_FE.elem_id()<< "    side_id"<<F_FE.side_id() <<std::endl;
         compute_face_dofs(mat,cont,QRuleFace::qp_weights,shape_face_.eval(), normal, F_FE);
         cont++;


    }


	   template<typename Mat1, typename Mat2>
	   void resize(Mat1& mat1,const Mat2& mat2)
	   {
	   	for(Integer i=0;i<NDofs;i++)
	   		mat1(i,i)=mat1(i,i)/mat2(i,i);
	   }


    template<typename Mat>
    void compute(Mat& mat, FiniteElem<Elem>& C_FE, FiniteElem<Elem>& F_FE)
    {
    	Integer cont=0;
    	// compute_aux<Cont>(mat,cont,C_FE,F_FE);
    	// compute_face_aux<0>(mat,cont,C_FE,F_FE);
    	// compute_momentum_aux<0>(mat,cont,C_FE,F_FE);
        // std::cout<<"compute first"<<std::endl;
    	for(Integer i=0;i<ManifoldDim+1;i++)
    		{F_FE.init_boundary(i,true);
    	     compute<entity[0]>(mat,cont,C_FE,F_FE);}
  
        cont=0;
        C_FE.init_inv_jac();
        F_FE.init(C_FE.elem_id());
     	for(Integer i=0;i<ManifoldDim+1;i++)
    		{F_FE.init_boundary(i,true);
    	     compute<entity[0]>(mat_coarse_,cont,C_FE,F_FE);}


    	// std::cout<<"pre mat"<<std::endl;
    	// std::cout<<mat<<std::endl;
    	// std::cout<<"mat_coarse_"<<std::endl;
    	// std::cout<<mat_coarse_<<std::endl;

        resize(mat,mat_coarse_);
    	// std::cout<<"mat"<<std::endl;
    	// std::cout<<mat<<std::endl;

    	// compute<entity[1]>(mat,C_FE,F_FE);
    }  


private:
	ShapeFace shape_face_;
	QPpointsFace qp_points_face_;
	QPpointsFace reference_qp_points_face_;
	Map C_map_;
	TraceMap trace_map_;
	Map F_map_;
	// Coefficients coeffs_;
	Array<Real,NDofs> coeff_;
	Matrix<Real,NDofs,NDofs> mat_coarse_;
};






template <Integer Dim,Integer ManifoldDim,Integer Continuity, Integer NComponents_,Integer QR>
class DofAux<ElementFunctionSpace<Simplex<Dim,ManifoldDim>,RaviartThomasFE,1,Continuity,NComponents_>,QR>
{
 public:

    using Elem=Simplex<Dim,ManifoldDim>;
    static constexpr Integer Order=1;
 	using Face=typename VolumeOrSurfaceElem<Elem,false>::type;
 	static constexpr Integer FEFamily=RaviartThomasFE;
 	static constexpr Integer NComponents=1; 
 	using RTnBaseFunctionSpace=BaseFunctionSpace<FEFamily,Order,Continuity,NComponents>;
 	using PnBaseFunctionSpace=BaseFunctionSpace<LagrangeFE,Order,Continuity,NComponents>;
 	using ElementFunctionSpace=ElementFunctionSpace<Elem,FEFamily,Order,Continuity,NComponents>;
 	static constexpr auto entity=ElementFunctionSpace::entity;
 	static constexpr auto dofs_per_entity=ElementFunctionSpace::dofs_per_entity;

    static constexpr Integer NDofs=FunctionSpaceDofsPerElem<ElementFunctionSpace>::value;

    static constexpr Integer NDofsFace= NComponents * dofs_per_entity[0] * ElemEntityCombinations<Elem,entity[0]>::value;
    static constexpr Integer NDofsMomentum= NComponents * dofs_per_entity[1] * ElemEntityCombinations<Elem,entity[1]>::value;

    // we have to integrate so + 1
    static constexpr Integer QRuleOrder=QuadratureOrder<IdentityOperator, RTnBaseFunctionSpace>::value+1;
    

    using Map=MapFromReference<IdentityOperator,Elem,FEFamily> ;
    using TraceMap=MapFromReference<TraceOperator,Face,FEFamily> ;

    using F_Lagrangian_Map=MapFromReference<IdentityOperator,Elem,LagrangeFE> ;


    using QRuleFaceReference=typename QuadratureRule<QR>:: template rule<Face,QRuleOrder>;
    using QRuleFace=Boundary2VolumetricQuadratureRule<QRuleFaceReference>;
    using ShapeFace=ShapeFunction<Elem,RTnBaseFunctionSpace,IdentityOperator,QRuleFace>;
    using ShapeFaceValue=typename ShapeFace::type;
    using QPpointsFace=typename ShapeFace::qp_points_type;

    // using ShapeTraceFace=ShapeFunction<Elem,RTnBaseFunctionSpace,TraceOperator,QRuleFaceReference>;

    using ShapeTraceFace=ShapeFunction<Elem,RTnBaseFunctionSpace,TraceOperator,QRuleFaceReference>;
    using trace_type=typename ShapeTraceFace::trace_type;

    using QRuleMomentum=typename QuadratureRule<QR>:: template rule<Elem,QRuleOrder>;
    using ShapeMomentum=ShapeFunction<Elem,RTnBaseFunctionSpace,IdentityOperator,QRuleMomentum>;
    using ShapeMomentumValue=typename ShapeMomentum::type;
    using QPpointsMomentum=typename ShapeMomentum::qp_points_type;

    using QRuleSimplexPointsReference=typename QuadratureRule<QR>:: template rule<Elem,2>;
    using ShapeSimplexPoints=ShapeFunction<Elem,RTnBaseFunctionSpace,IdentityOperator,QRuleSimplexPointsReference>;



	using QRuleBarycenterReference=typename QuadratureRule<QR>:: template rule<Elem,1>;
    using ShapeBarycenter=ShapeFunction<Elem,RTnBaseFunctionSpace,IdentityOperator,QRuleBarycenterReference>;
    using QPpointsBarycenter=typename ShapeBarycenter::qp_points_type;


    using Coefficients=ShapeFunctionCoefficientsCollectionSingleType<Elem,FEFamily>;
    using ShapeLagrangianMomentum=ShapeFunction<Elem,PnBaseFunctionSpace,IdentityOperator,QRuleMomentum>;

    template<typename Output, typename Weight, typename ShapeValue,typename ShapeTraceValue,typename Normal>
    void compute_face_dofs(Output& out, Integer& fine_dof, const Weight& weights,const ShapeValue& shape_value,ShapeTraceValue& shape_trace_value, const Normal& normal, FiniteElem<Elem>& FE)
    {
     const auto& NQPoints=ShapeValue::NQPoints;
     const auto& Ndofs=Output::Rows;

     // std::cout<<"compute_face_dofs"<<std::endl;


     // std::cout<<"normal"<<std::endl;
     // std::cout<<normal<<std::endl;
   
     for(Integer coarse_dof=0;coarse_dof<NDofs;coarse_dof++)
     {
        out(fine_dof,coarse_dof)=0;

     	for(Integer qp=0;qp<NQPoints;qp++)
     	{
     		for(Integer i=0;i<ManifoldDim;i++)
     		{
     			out(fine_dof,coarse_dof)+= weights[qp]*shape_value[coarse_dof][qp][i]*normal[i]*shape_trace_value[qp](0,0);
     		}
     	}
       out(fine_dof,coarse_dof)*=FE.get_det_side();

     }
     // std::cout<<"out"<<std::endl;
     // std::cout<<out<<std::endl;    
    }

    // template<typename Output, typename Weight, typename ShapeValue,typename Normal>
    // void compute_face_dofs(Output& out, Integer& fine_dof, const Weight& weights,const ShapeValue& shape_value, const Normal& normal, FiniteElem<Elem>& FE)
    // {
    //  const auto& NQPoints=ShapeValue::NQPoints;
    //  const auto& Ndofs=Output::Rows;

    //  // std::cout<<"compute_face_dofs"<<std::endl;
    //  // std::cout<<"shape_value"<<std::endl;
    //  // std::cout<<shape_value<<std::endl; 

    //  // std::cout<<"normal"<<std::endl;
    //  // std::cout<<normal<<std::endl;
    //  // std::cout<<"weights"<<std::endl;
    //  // std::cout<<weights<<std::endl;
    //  // std::cout<<"FE.get_det_side()"<<std::endl;
    //  // std::cout<<FE.get_det_side()<<std::endl;
    //  for(Integer coarse_dof=0;coarse_dof<NDofs;coarse_dof++)
    //  {
    //     out(fine_dof,coarse_dof)=0;

    //  	for(Integer qp=0;qp<NQPoints;qp++)
    //  	{
    //  		for(Integer i=0;i<ManifoldDim;i++)
    //  		{
    //  			out(fine_dof,coarse_dof)+= weights[qp]*shape_value[coarse_dof][qp][i]*normal[i];
    //  		}
    //  	}
    //    out(fine_dof,coarse_dof)*=FE.get_det_side();

    //  }
    //  // std::cout<<"out"<<std::endl;
    //  // std::cout<<out<<std::endl;    
    // }



    // template<typename Output, typename Weight, typename ShapeValue, typename ShapeLagrangeValue>
    // void compute_momentum_dofs(Output& out, Integer& fine_dof,  const Weight& weights, ShapeValue& shape_value,const Integer component, ShapeLagrangeValue& shape_lagrange_value,const FiniteElem<Elem>& FE)
    // {
    //  const auto& NQPoints=ShapeValue::NQPoints;
    //  const auto& Ndofs=Output::Rows;
     

    //   // std::cout<<"compute_momentum_dofs"<<std::endl;
    //   // std::cout<<"weights"<<std::endl;
    //   // std::cout<<weights<<std::endl;

    //   // std::cout<<"shape_value"<<std::endl;
    //   // std::cout<<shape_value<<std::endl;

    //   // std::cout<<"shape_lagrange_value"<<std::endl;
    //   // std::cout<<shape_lagrange_value<<std::endl;

    //   // std::cout<<"FE.get_det()"<<std::endl;
    //   // std::cout<<FE.get_det()<<std::endl;

    //  for(Integer coarse_dof=0;coarse_dof<Ndofs;coarse_dof++)
    //  {
    //     out(fine_dof,coarse_dof)=0;
    //     // std::cout<<"n_dof="<<n_dof<<std::endl;
    //  	for(Integer qp=0;qp<NQPoints;qp++)
    //  	{
    //  			out(fine_dof,coarse_dof)+= weights[qp]*shape_value[coarse_dof][qp][component]*shape_lagrange_value[qp](0,0);
    //  			// std::cout<<out(n_dof,cont)<<" ";
    //  	}

    //    // std::cout<<out(n_dof,cont)<<std::endl;
    //    out(fine_dof,coarse_dof)*=FE.get_det();

    //  }
    //  // std::cout<<"out"<<std::endl;
    //  // std::cout<<out<<std::endl;    
    // }



    template<typename Output, typename ShapeValue>
    void compute_barycenter_dofs(Output& out, Integer& fine_dof,  const Integer component, ShapeValue& shape_value)
    {
     const auto& NQPoints=ShapeValue::NQPoints;
     const auto& Ndofs=Output::Rows;


     // std::cout<<"Ndofs="<<Ndofs<<" NQPoints="<<NQPoints<<std::endl;


     for(Integer coarse_dof=0;coarse_dof<Ndofs;coarse_dof++)
     {
        out(fine_dof,coarse_dof)=0;
        // std::cout<<"fine_dof="<<fine_dof<<" coarse_dof="<<coarse_dof<<std::endl;
     	for(Integer qp=0;qp<1;qp++)
     	{
     			out(fine_dof,coarse_dof)=shape_value[coarse_dof][qp][component];
     			// std::cout<<out(fine_dof,coarse_dof)<<" ";
     	}
     	// std::cout<<std::endl;

       // std::cout<<out(n_dof,cont)<<std::endl;

     }
  
    }


    template<typename Output, typename ShapeValue,typename Normal>
    void compute_simplex_nodes_dofs(Output& out, Integer& fine_dof,const Integer& k, const Integer& local_node, const Normal& normal, const ShapeValue& shape_value)
    {
     const auto& NQPoints=ShapeValue::NQPoints;
     const auto& Ndofs=Output::Rows;

     // std::cout<<"compute_simplex_nodes_dofs="<<std::endl;
     // std::cout<<"Ndofs="<<Ndofs<<" NQPoints="<<NQPoints<<std::endl;
     // std::cout<<"fine_dof="<<fine_dof<<std::endl;
     // std::cout<<"k="<<k<<std::endl;
     // std::cout<<"local_node="<<local_node<<std::endl;

     // Integer coarse_dof=k*ManifoldDim+local_node;
     for(Integer coarse_dof=0;coarse_dof<Ndofs;coarse_dof++)
     {
        out(fine_dof,coarse_dof)=0;
        // std::cout<<" coarse_dof="<<coarse_dof<<std::endl;
        // std::cout<<shape_value[coarse_dof]<<std::endl;
     	for(Integer comp=0;comp<ManifoldDim;comp++)
     	{
     			out(fine_dof,coarse_dof)+=shape_value[coarse_dof][local_node][comp] * normal[comp];
     			
     	}
     	// std::cout<<out(fine_dof,coarse_dof)<<" ";
     	// std::cout<<std::endl;

       // std::cout<<out(n_dof,cont)<<std::endl;

     }


    }



    // // face dofs
    // template<Integer EntityDim, typename Mat>
    // std::enable_if_t< (EntityDim==entity[0]), void>
    // compute(Mat& mat, Integer& cont, const FiniteElem<Elem>& C_FE, FiniteElem<Elem>& F_FE)
    // {



    //  	auto mesh_ptr=F_FE.mesh_ptr();
    // 	const auto& signed_normal=mesh_ptr->signed_normal();
    //     C_map_.init(C_FE);
    //     shape_face_.init_map(C_map_);
    //     // coeffs_.init(*mesh_ptr);

    //     SingleShapeFunctionCoefficientsCollection<Elem,FEFamily,Order>::apply(signed_normal.sign(C_FE.elem_id()),coeff_);
    // 	F_FE.side_transform_point(qp_points_face_,QRuleFaceReference::qp_points);
    //  // 	std::cout<<"QRuleFaceReference::qp_points"<<std::endl;
    // 	// std::cout<<QRuleFaceReference::qp_points<<std::endl;
    // 	// std::cout<<"qp_points_face_"<<std::endl;
    // 	// std::cout<<qp_points_face_<<std::endl;
    //  // transformed_point_row_=C_FE_.inv_jac()*transformed_point_col_;
    //     C_FE.points_to_reference_points(reference_qp_points_face_,qp_points_face_);
    //     // std::cout<<"qui7 "<<std::endl;
    // 	// std::cout<<"reference_qp_points_face_"<<std::endl;
    // 	std::cout<<reference_qp_points_face_<<std::endl;
    //     shape_face_.init(coeff_,reference_qp_points_face_,C_FE);
    //     // // Vector<Real,NDofsFace> out; 

    //     const auto& normal=signed_normal.normals()[F_FE.elem_id()][F_FE.side_id()];

    //     Integer id=F_FE.side_id();
    //     // std::cout<<"qui2"<<std::endl;
    //     trace_map_.init(F_FE);
    //     // std::cout<<"qui3"<<std::endl;
    //     shape_trace_face_.init_map(trace_map_);
    //     // std::cout<<"qui4"<<std::endl;
    //     shape_trace_face_.init(coeff_,id,C_FE);

        
    //     std::cout<<"reference_qp_points_face_"<<std::endl;
    //     std::cout<<reference_qp_points_face_<<std::endl;
    //     std::cout<<"normal"<<std::endl;
    //     std::cout<<normal<<std::endl;
    //     std::cout<<"shape_face_.eval()"<<std::endl;
    //     std::cout<<shape_face_.eval()<<std::endl;
    //     std::cout<<"shape_trace_face_.eval()"<<std::endl;
    //     std::cout<<shape_trace_face_.eval()<<std::endl;

    //     for(Integer i=0;i<dofs_per_entity[0]; i++)
    //     {
    //      compute_face_dofs(mat,cont,QRuleFace::qp_weights,shape_face_.eval() ,shape_trace_face_.eval()[i], normal, F_FE);
    //      cont++;
    //     }





    //  // F_map_.init(F_FE_);
    // 	// QPpointsFace frfrfr(5);
    // 	// typename QRuleFaceReference::qp_points_type feef(5);
    // 	// // std::cout<<"compute "<<std::endl;
    // 	// auto mesh_ptr=F_FE.mesh_ptr();
    // 	// // std::cout<<"qui1 "<<std::endl;
    // 	// const auto& signed_normal=mesh_ptr->signed_normal();
    // 	// // std::cout<<"qui2 "<<std::endl;
    //  //    C_map_.init(C_FE);
    //  //    // std::cout<<"qui3 "<<std::endl;
    //  //    shape_face_.init_map(C_map_);
    //  //    // std::cout<<"qui4 "<<std::endl;
    //  //    // coeffs_.init(*mesh_ptr);

    //  //    SingleShapeFunctionCoefficientsCollection<Elem,FEFamily,Order>::apply(signed_normal.sign(C_FE.elem_id()),coeff_);
    //  //    // std::cout<<"qui5 "<<std::endl;
    // 	// F_FE.side_transform_point(qp_points_face_,QRuleFaceReference::qp_points);
    // 	// // std::cout<<"qui6 "<<std::endl;
    //  // // 	std::cout<<"QRuleFaceReference::qp_points"<<std::endl;
    // 	// // std::cout<<QRuleFaceReference::qp_points<<std::endl;
    // 	// // std::cout<<"qp_points_face_"<<std::endl;
    // 	// // std::cout<<qp_points_face_<<std::endl;
    //  // // transformed_point_row_=C_FE_.inv_jac()*transformed_point_col_;
    //  //    C_FE.points_to_reference_points(reference_qp_points_face_,qp_points_face_);
    //  //    // std::cout<<"qui7 "<<std::endl;
    // 	// // std::cout<<"reference_qp_points_face_"<<std::endl;
    // 	// // std::cout<<reference_qp_points_face_<<std::endl;
    //  //    shape_face_.init(coeff_,reference_qp_points_face_);
    //  //    std::cout<<"shape_face_.eval()"<<std::endl;
    //  //    // std::cout<<shape_face_.eval()<<std::endl;
    //  //    // Vector<Real,NDofsFace> out; 

    //  //    const auto& normal=signed_normal.normals()[F_FE.elem_id()][F_FE.side_id()];
        
    //  //    // std::cout<<"qui1"<<std::endl;
    //  //    Integer id=F_FE.side_id();
    //  //    // std::cout<<"qui2"<<std::endl;
    //  //    trace_map_.init(F_FE);
    //  //    // std::cout<<"qui3"<<std::endl;
    //  //    shape_trace_face_.init_map(trace_map_);
    //  //    // std::cout<<"qui4"<<std::endl;
    //  //    shape_trace_face_.init(coeff_,id);
    //  //    // std::cout<<"F_FE.side_id()"<<std::endl;
    //  //    // std::cout<<F_FE.side_id()<<std::endl;

    //  //    // std::cout<<"shape_trace_face_.eval()"<<std::endl;
    //  //    // std::cout<<shape_trace_face_.eval()<<std::endl;
        
    //  //    for(Integer i=0;i<dofs_per_entity[0]; i++)
    //  //    {
    //  //     compute_face_dofs(mat,cont,QRuleFace::qp_weights,shape_face_.eval() ,shape_trace_face_.eval()[i], normal, F_FE);
    //  //     cont++;
    //  //    }
    // 	std::cout<<"mat"<<std::endl;
    // 	std::cout<<mat<<std::endl;

 


    //  // // qp_points_face_
    //  // // shape_face_.init_map(C_map);
    // }


    // // momentum dofs
    // template<Integer EntityDim, typename Mat>
    // std::enable_if_t< (EntityDim==entity[1]), void> 
    // compute(Mat& mat, Integer& cont, const FiniteElem<Elem>& C_FE, FiniteElem<Elem>& F_FE)
    // {
    //  // qp_points_momentum_
    //  auto mesh_ptr=F_FE.mesh_ptr();
    //  const auto& signed_normal=mesh_ptr->signed_normal();
    //  C_map_.init(C_FE);
    //  F_Pn_map_.init(F_FE);
    //  shape_momentum_.init_map(C_map_);
    //  shape_lagrange_momentum_.init_map(F_Pn_map_);
    //  SingleShapeFunctionCoefficientsCollection<Elem,FEFamily,Order>::apply(signed_normal.sign(C_FE.elem_id()),coeff_);
    //  F_FE.transform_point(qp_points_momentum_,QRuleMomentum::qp_points);
    //  C_FE.points_to_reference_points(reference_qp_points_momentum_,qp_points_momentum_);

    //  std::cout<<"QRuleMomentum::qp_points"<<std::endl;
    //  // std::cout<<QRuleMomentum::qp_points<<std::endl;
    //  // std::cout<<"qp_points_momentum_"<<std::endl;
    //  // std::cout<<qp_points_momentum_<<std::endl;
    //  // std::cout<<"reference_qp_points_momentum_"<<std::endl;
    //  // std::cout<<reference_qp_points_momentum_<<std::endl;

    //  shape_momentum_.init(reference_qp_points_momentum_);
    //  shape_lagrange_momentum_.init(reference_qp_points_momentum_);
     
    //  int cont2=0;
    //  // Vector<Real,NDofsMomentum > out; 
    //  for(Integer i=NDofsFace;i<NDofs;i++)
    //  {
    //  compute_momentum_dofs(mat,cont,QRuleFace::qp_weights,shape_momentum_.eval(),cont2,shape_lagrange_momentum_.eval()[cont2], F_FE);
    //  cont++;
    //  cont2++;

    //  }


    //  // std::cout<<"QRuleMomentum::qp_points"<<std::endl;
    //  // std::cout<<QRuleMomentum::qp_points<<std::endl;

    //  // std::cout<<"qp_points_momentum_"<<std::endl;
    //  // std::cout<<qp_points_momentum_<<std::endl;



    //  // std::cout<<"shape_momentum_.eval()"<<std::endl;
    //  // std::cout<<shape_momentum_.eval()<<std::endl;
    //  // std::cout<<"shape_lagrange_momentum_.eval()"<<std::endl;
    //  // std::cout<<shape_lagrange_momentum_.eval()<<std::endl;

    //  // std::cout<<"mat"<<std::endl;
    //  // std::cout<<mat<<std::endl;
    // }



    static constexpr const Array<Array<Real, ManifoldDim>,ManifoldDim+1> reference_nodes_arr=Elem::reference();
    static constexpr const auto reference_nodes=ArrayOfArray2Matrix(reference_nodes_arr);


	    template<Integer EntityDim, typename Mat>
	    std::enable_if_t< (EntityDim==entity[0]), void> 
	    compute(Mat& mat, Integer& cont, FiniteElem<Elem>& C_FE, FiniteElem<Elem>& F_FE)
	    {





     	auto mesh_ptr=F_FE.mesh_ptr();
    	const auto& signed_normal=mesh_ptr->signed_normal();
     
    	F_FE.init(F_FE.elem_id());
    	C_FE.init(C_FE.elem_id());
    	C_FE.init_inv_jac();
        C_map_.init(C_FE);
        // shape_simplex_points_.init_map(C_map_);


  
        FQPValues<Matrix<Real,ManifoldDim,1>,ManifoldDim+1,NDofs> func_values;


        SingleShapeFunctionCoefficientsCollection<Elem,FEFamily,Order>::apply(signed_normal.sign(C_FE.elem_id()),coeff_);
     	F_FE.transform_point(qp_points_simplex_,reference_nodes);
        C_FE.points_to_reference_points(reference_qp_points_simplex_,qp_points_simplex_);


        DynamicShapeFunctionValue<ElementFunctionSpace,IdentityOperator>::apply4(func_values,reference_qp_points_simplex_,C_FE,C_map_(),coeff_); 




        // coeffs_.init(*mesh_ptr);

        // for(Integer i=0;i<ManifoldDim+1;i++)
    	// F_FE.transform_point(qp_points_simplex_,reference_nodes);
     // 	std::cout<<"QRuleFaceReference::qp_points"<<std::endl;
    	// std::cout<<QRuleFaceReference::qp_points<<std::endl;
    	// std::cout<<"qp_points_face_"<<std::endl;
    	// std::cout<<qp_points_face_<<std::endl;
     // transformed_point_row_=C_FE_.inv_jac()*transformed_point_col_;
        // C_FE.points_to_reference_points(reference_qp_points_simplex_,qp_points_simplex_);
        // std::cout<<"qui7 "<<std::endl;
        std::cout<<"reference_nodes"<<std::endl;
        std::cout<<reference_nodes<<std::endl;
        std::cout<<"qp_points_simplex_"<<std::endl;
        std::cout<<qp_points_simplex_<<std::endl;
    	std::cout<<"reference_qp_points_simplex_"<<std::endl;
    	std::cout<<reference_qp_points_simplex_<<std::endl;
        // shape_simplex_points_.init(coeff_,reference_qp_points_simplex_,C_FE);
        // // Vector<Real,NDofsFace> out; 
    	// std::cout<<"interp func_values"<<std::endl;
    	// std::cout<<func_values<<std::endl;
        const auto& normal=signed_normal.normals()[F_FE.elem_id()][F_FE.side_id()];
        Integer local_nodes[ManifoldDim];


        Combinations<ManifoldDim + 1,ManifoldDim>::generate(F_FE.side_id(),local_nodes);


        for(Integer k=0;k<ManifoldDim;k++)
        {       	
     //    	std::cout<<"local_nodes"<<std::endl;
     //    	for(Integer m=0;m<ManifoldDim;m++)
     //    	std::cout<<local_nodes[m]<<" ";
     //    std::cout<<"local_nodes[k]="<<local_nodes[k]<<std::endl;
     //    std::cout<<std::endl;
    	// std::cout<<"func_values"<<std::endl;
    	// std::cout<<func_values<<std::endl; 
    	// std::cout<<"normal"<<std::endl;
    	// std::cout<<normal<<std::endl;        
         compute_simplex_nodes_dofs(mat,cont,k,local_nodes[k],normal,func_values);
         cont++;
        }


    	// std::cout<<"mat"<<std::endl;
    	// std::cout<<mat<<std::endl;
    }













	    template<Integer EntityDim, typename Mat>
	    std::enable_if_t< (EntityDim==entity[1]), void> 
	    compute(Mat& mat, Integer& cont, FiniteElem<Elem>& C_FE, FiniteElem<Elem>& F_FE)
	    {





     	auto mesh_ptr=F_FE.mesh_ptr();
    	const auto& signed_normal=mesh_ptr->signed_normal();
     
    	F_FE.init(F_FE.elem_id());
    	C_FE.init(C_FE.elem_id());
    	C_FE.init_inv_jac();
        C_map_.init(C_FE);
        // shape_simplex_points_.init_map(C_map_);


  
        FQPValues<Matrix<Real,ManifoldDim,1>,ManifoldDim+1,NDofs> func_values;


        SingleShapeFunctionCoefficientsCollection<Elem,FEFamily,Order>::apply(signed_normal.sign(C_FE.elem_id()),coeff_);
    	F_FE.transform_point(qp_points_barycenter_,QRuleBarycenterReference::qp_points);
         C_FE.points_to_reference_points(reference_qp_points_barycenter_,qp_points_barycenter_);


        DynamicShapeFunctionValue<ElementFunctionSpace,IdentityOperator>::apply4(func_values,reference_qp_points_barycenter_,C_FE,C_map_(),coeff_); 




        // coeffs_.init(*mesh_ptr);

        // for(Integer i=0;i<ManifoldDim+1;i++)
    	// F_FE.transform_point(qp_points_simplex_,reference_nodes);
     // 	std::cout<<"QRuleFaceReference::qp_points"<<std::endl;
    	// std::cout<<QRuleFaceReference::qp_points<<std::endl;
    	// std::cout<<"qp_points_face_"<<std::endl;
    	// std::cout<<qp_points_face_<<std::endl;
     // transformed_point_row_=C_FE_.inv_jac()*transformed_point_col_;
        // C_FE.points_to_reference_points(reference_qp_points_simplex_,qp_points_simplex_);
        // std::cout<<"qui7 "<<std::endl;
        std::cout<<"QRuleBarycenterReference::qp_points"<<std::endl;
        std::cout<<QRuleBarycenterReference::qp_points<<std::endl;
        std::cout<<"qp_points_barycenter_"<<std::endl;
        std::cout<<qp_points_barycenter_<<std::endl;
    	std::cout<<"reference_qp_points_barycenter_"<<std::endl;
    	std::cout<<reference_qp_points_barycenter_<<std::endl;
        // shape_simplex_points_.init(coeff_,reference_qp_points_simplex_,C_FE);
     //    // // Vector<Real,NDofsFace> out; 
    	// std::cout<<"interp func_values"<<std::endl;
    	// std::cout<<func_values<<std::endl;
        // const auto& normal=signed_normal.normals()[F_FE.elem_id()][F_FE.side_id()];
        // Integer local_nodes[ManifoldDim];


        // Combinations<ManifoldDim + 1,ManifoldDim>::generate(F_FE.side_id(),local_nodes);


        for(Integer k=0;k<dofs_per_entity[1];k++)
        {
         compute_barycenter_dofs(mat,cont,k,func_values);
         cont++;
        }


    	// std::cout<<"mat"<<std::endl;
    	// std::cout<<mat<<std::endl;
    }



    // template<Integer EntityDim, typename Mat>
    // std::enable_if_t< (EntityDim==entity[1]), void> 
    // compute(Mat& mat, Integer& cont, const FiniteElem<Elem>& C_FE, FiniteElem<Elem>& F_FE)
    // {




    //  	auto mesh_ptr=F_FE.mesh_ptr();
    // 	const auto& signed_normal=mesh_ptr->signed_normal();
    //     C_map_.init(C_FE);
    //     shape_barycenter_.init_map(C_map_);
    //     // coeffs_.init(*mesh_ptr);

    //     SingleShapeFunctionCoefficientsCollection<Elem,FEFamily,Order>::apply(signed_normal.sign(C_FE.elem_id()),coeff_);
    // 	F_FE.transform_point(qp_points_barycenter_,QRuleBarycenterReference::qp_points);
    //  // 	std::cout<<"QRuleFaceReference::qp_points"<<std::endl;
    // 	// std::cout<<QRuleFaceReference::qp_points<<std::endl;
    // 	// std::cout<<"qp_points_face_"<<std::endl;
    // 	// std::cout<<qp_points_face_<<std::endl;
    //  // transformed_point_row_=C_FE_.inv_jac()*transformed_point_col_;
    //     C_FE.points_to_reference_points(reference_qp_points_barycenter_,qp_points_barycenter_);
    //     // std::cout<<"qui7 "<<std::endl;
    //  //    std::cout<<"QRuleBarycenterReference::qp_points"<<std::endl;
    //  //    std::cout<<QRuleBarycenterReference::qp_points<<std::endl;
    //  //    std::cout<<"qp_points_barycenter_"<<std::endl;
    //  //    std::cout<<qp_points_barycenter_<<std::endl;
    // 	// std::cout<<"reference_qp_points_barycenter_"<<std::endl;
    // 	// std::cout<<reference_qp_points_barycenter_<<std::endl;
    //     shape_barycenter_.init(coeff_,reference_qp_points_barycenter_,C_FE);
    //     // // Vector<Real,NDofsFace> out; 

    //     const auto& normal=signed_normal.normals()[F_FE.elem_id()][F_FE.side_id()];
    //     for(Integer k=0;k<dofs_per_entity[1];k++)
    //     {
    //      compute_barycenter_dofs(mat,cont,k,shape_barycenter_.eval());
    //      cont++;
    //     }





    //  // auto mesh_ptr=F_FE.mesh_ptr();
    //  // const auto& signed_normal=mesh_ptr->signed_normal();
    //  // C_map_.init(C_FE);
    //  // F_Pn_map_.init(F_FE);
    //  // shape_momentum_.init_map(C_map_);
    //  // shape_lagrange_momentum_.init_map(F_Pn_map_);
    //  // SingleShapeFunctionCoefficientsCollection<Elem,FEFamily,Order>::apply(signed_normal.sign(C_FE.elem_id()),coeff_);
    //  // F_FE.transform_point(qp_points_momentum_,QRuleMomentum::qp_points);
    //  // C_FE.points_to_reference_points(reference_qp_points_momentum_,qp_points_momentum_);

    //  // std::cout<<"QRuleMomentum::qp_points"<<std::endl;
    //  // // std::cout<<QRuleMomentum::qp_points<<std::endl;
    //  // // std::cout<<"qp_points_momentum_"<<std::endl;
    //  // // std::cout<<qp_points_momentum_<<std::endl;
    //  // // std::cout<<"reference_qp_points_momentum_"<<std::endl;
    //  // // std::cout<<reference_qp_points_momentum_<<std::endl;

    //  // shape_momentum_.init(reference_qp_points_momentum_);
    //  // shape_lagrange_momentum_.init(reference_qp_points_momentum_);
     
    //  // int cont2=0;
    //  // // Vector<Real,NDofsMomentum > out; 
    //  // for(Integer i=NDofsFace;i<NDofs;i++)
    //  // {
    //  // compute_momentum_dofs(mat,cont,QRuleFace::qp_weights,shape_momentum_.eval(),cont2,shape_lagrange_momentum_.eval()[cont2], F_FE);
    //  // cont++;
    //  // cont2++;

    //  // }

    // }

    template<typename Mat>
    void compute(Mat& mat, FiniteElem<Elem>& C_FE, FiniteElem<Elem>& F_FE)
    {
    	Integer cont=0;
    	// compute_aux<Cont>(mat,cont,C_FE,F_FE);
    	// compute_face_aux<0>(mat,cont,C_FE,F_FE);
    	// compute_momentum_aux<0>(mat,cont,C_FE,F_FE);
        std::cout<<"compute first"<<std::endl;

        // C_FE_new.init(C_FE.elem_id());

        

       
        
    	for(Integer i=0;i<ManifoldDim+1;i++)
    		{
    			F_FE.init_boundary(i,true);
    			compute<entity[0]>(mat,cont,C_FE,F_FE);
    	 }

    	 std::cout<<"compute second"<<std::endl;

    	compute<entity[1]>(mat,cont,C_FE,F_FE);

    	


     std::cout<<"last final mat"<<std::endl;
     std::cout<<mat<<std::endl;

    	// compute<entity[1]>(mat,C_FE,F_FE);
    }  


private:
	ShapeBarycenter shape_barycenter_;
	ShapeFace shape_face_;
	ShapeTraceFace shape_trace_face_;
	ShapeMomentum shape_momentum_;
	ShapeLagrangianMomentum shape_lagrange_momentum_;

	ShapeSimplexPoints shape_simplex_points_;

	QPpointsBarycenter qp_points_barycenter_;
	QPpointsBarycenter reference_qp_points_barycenter_;
    
	Matrix<Real,ManifoldDim+1,ManifoldDim> qp_points_simplex_;
	Matrix<Real,ManifoldDim+1,ManifoldDim> reference_qp_points_simplex_;

	QPpointsFace qp_points_face_;
	QPpointsFace reference_qp_points_face_;
	QPpointsMomentum qp_points_momentum_;
	QPpointsMomentum reference_qp_points_momentum_;
	Map C_map_;
	TraceMap trace_map_;
	Map F_map_;
	F_Lagrangian_Map F_Pn_map_;
	// Coefficients coeffs_;
	Array<Real,NDofs> coeff_;
};



// template <Integer Dim,Integer ManifoldDim,Integer Continuity, Integer NComponents_,Integer QR>
// class DofAux<ElementFunctionSpace<Simplex<Dim,ManifoldDim>,RaviartThomasFE,1,Continuity,NComponents_>,QR>
// {
//  public:

//     using Elem=Simplex<Dim,ManifoldDim>;
//     static constexpr Integer Order=1;
//  	using Face=typename VolumeOrSurfaceElem<Elem,false>::type;
//  	static constexpr Integer FEFamily=RaviartThomasFE;
//  	static constexpr Integer NComponents=1; 
//  	using RTnBaseFunctionSpace=BaseFunctionSpace<FEFamily,Order,Continuity,NComponents>;
//  	using PnBaseFunctionSpace=BaseFunctionSpace<LagrangeFE,Order,Continuity,NComponents>;
//  	using ElementFunctionSpace=ElementFunctionSpace<Elem,FEFamily,Order,Continuity,NComponents>;
//  	static constexpr auto entity=ElementFunctionSpace::entity;
//  	static constexpr auto dofs_per_entity=ElementFunctionSpace::dofs_per_entity;

//     static constexpr Integer NDofs=FunctionSpaceDofsPerElem<ElementFunctionSpace>::value;

//     static constexpr Integer NDofsFace= NComponents * dofs_per_entity[0] * ElemEntityCombinations<Elem,entity[0]>::value;
//     static constexpr Integer NDofsMomentum= NComponents * dofs_per_entity[1] * ElemEntityCombinations<Elem,entity[1]>::value;

//     // we have to integrate so + 1
//     static constexpr Integer QRuleOrder=QuadratureOrder<IdentityOperator, RTnBaseFunctionSpace>::value+1;
    

//     using Map=MapFromReference<IdentityOperator,Elem,FEFamily> ;
//     using TraceMap=MapFromReference<TraceOperator,Face,FEFamily> ;

//     using F_Lagrangian_Map=MapFromReference<IdentityOperator,Elem,LagrangeFE> ;


//     using QRuleFaceReference=typename QuadratureRule<QR>:: template rule<Face,QRuleOrder>;
//     using QRuleFace=Boundary2VolumetricQuadratureRule<QRuleFaceReference>;
//     using ShapeFace=ShapeFunction<Elem,RTnBaseFunctionSpace,IdentityOperator,QRuleFace>;
//     using ShapeFaceValue=typename ShapeFace::type;
//     using QPpointsFace=typename ShapeFace::qp_points_type;


//     using ShapeTraceFace=ShapeFunction<Elem,RTnBaseFunctionSpace,TraceOperator,QRuleFaceReference>;
//     using trace_type=typename ShapeTraceFace::trace_type;

//     using QRuleMomentum=typename QuadratureRule<QR>:: template rule<Elem,QRuleOrder>;
//     using ShapeMomentum=ShapeFunction<Elem,RTnBaseFunctionSpace,IdentityOperator,QRuleMomentum>;
//     using ShapeMomentumValue=typename ShapeMomentum::type;
//     using QPpointsMomentum=typename ShapeMomentum::qp_points_type;

//     using Coefficients=ShapeFunctionCoefficientsCollectionSingleType<Elem,FEFamily>;
//     using ShapeLagrangianMomentum=ShapeFunction<Elem,PnBaseFunctionSpace,IdentityOperator,QRuleMomentum>;

//   // template<typename Mat,typename Values>
//   // static void compute(Mat& local_mat, Values& func_values )
//   // {
//   // 	local_mat(M,N)=func_values[N][0](0,0);
//   // }
//     template<typename Output, typename Weight, typename ShapeValue,typename ShapeTraceValue,typename Normal>
//     void compute_face_dofs(Output& out, Integer& fine_dof, const Weight& weights,const ShapeValue& shape_value,ShapeTraceValue& shape_trace_value, const Normal& normal, FiniteElem<Elem>& FE)
//     {
//      const auto& NQPoints=ShapeValue::NQPoints;
//      const auto& Ndofs=Output::Rows;

//      // std::cout<<"compute_face_dofs"<<std::endl;


//      // std::cout<<"normal"<<std::endl;
//      // std::cout<<normal<<std::endl;
   
//      for(Integer coarse_dof=0;coarse_dof<NDofs;coarse_dof++)
//      {
//         out(fine_dof,coarse_dof)=0;

//      	for(Integer qp=0;qp<NQPoints;qp++)
//      	{
//      		for(Integer i=0;i<ManifoldDim;i++)
//      		{
//      			out(fine_dof,coarse_dof)+= weights[qp]*shape_value[coarse_dof][qp][i]*normal[i]*shape_trace_value[qp](0,0);
//      		}
//      	}
//        out(fine_dof,coarse_dof)*=FE.get_det_side();

//      }
//      // std::cout<<"out"<<std::endl;
//      // std::cout<<out<<std::endl;    
//     }

//     template<typename Output, typename Weight, typename ShapeValue, typename ShapeLagrangeValue>
//     void compute_momentum_dofs(Output& out, Integer& fine_dof,  const Weight& weights, ShapeValue& shape_value,const Integer component, ShapeLagrangeValue& shape_lagrange_value,const FiniteElem<Elem>& FE)
//     {
//      const auto& NQPoints=ShapeValue::NQPoints;
//      const auto& Ndofs=Output::Rows;
     

//       // std::cout<<"compute_momentum_dofs"<<std::endl;
//       // std::cout<<"weights"<<std::endl;
//       // std::cout<<weights<<std::endl;

//       // std::cout<<"shape_value"<<std::endl;
//       // std::cout<<shape_value<<std::endl;

//       // std::cout<<"shape_lagrange_value"<<std::endl;
//       // std::cout<<shape_lagrange_value<<std::endl;

//       // std::cout<<"FE.get_det()"<<std::endl;
//       // std::cout<<FE.get_det()<<std::endl;

//      for(Integer coarse_dof=0;coarse_dof<Ndofs;coarse_dof++)
//      {
//         out(fine_dof,coarse_dof)=0;
//         // std::cout<<"n_dof="<<n_dof<<std::endl;
//      	for(Integer qp=0;qp<NQPoints;qp++)
//      	{
//      			out(fine_dof,coarse_dof)+= weights[qp]*shape_value[coarse_dof][qp][component]*shape_lagrange_value[qp](0,0);
//      			// std::cout<<out(n_dof,cont)<<" ";
//      	}

//        // std::cout<<out(n_dof,cont)<<std::endl;
//        out(fine_dof,coarse_dof)*=FE.get_det();

//      }
//      // std::cout<<"out"<<std::endl;
//      // std::cout<<out<<std::endl;    
//     }

//     // face dofs
//     template<Integer EntityDim, typename Mat>
//     std::enable_if_t< (EntityDim==entity[0]), void>
//     compute(Mat& mat, Integer& cont, const FiniteElem<Elem>& C_FE, FiniteElem<Elem>& F_FE)
//     {
//      // F_map_.init(F_FE_);
//     	// QPpointsFace frfrfr(5);
//     	// typename QRuleFaceReference::qp_points_type feef(5);
//     	// std::cout<<"compute "<<std::endl;
//     	auto mesh_ptr=F_FE.mesh_ptr();
//     	// std::cout<<"qui1 "<<std::endl;
//     	const auto& signed_normal=mesh_ptr->signed_normal();
//     	// std::cout<<"qui2 "<<std::endl;
//         C_map_.init(C_FE);
//         // std::cout<<"qui3 "<<std::endl;
//         shape_face_.init_map(C_map_);
//         // std::cout<<"qui4 "<<std::endl;
//         // coeffs_.init(*mesh_ptr);

//         SingleShapeFunctionCoefficientsCollection<Elem,FEFamily,Order>::apply(signed_normal.sign(C_FE.elem_id()),coeff_);
//         // std::cout<<"qui5 "<<std::endl;
//     	F_FE.side_transform_point(qp_points_face_,QRuleFaceReference::qp_points);
//     	// std::cout<<"qui6 "<<std::endl;
//      // 	std::cout<<"QRuleFaceReference::qp_points"<<std::endl;
//     	// std::cout<<QRuleFaceReference::qp_points<<std::endl;
//     	// std::cout<<"qp_points_face_"<<std::endl;
//     	// std::cout<<qp_points_face_<<std::endl;
//      // transformed_point_row_=C_FE_.inv_jac()*transformed_point_col_;
//         C_FE.points_to_reference_points(reference_qp_points_face_,qp_points_face_);
//         // std::cout<<"qui7 "<<std::endl;
//     	// std::cout<<"reference_qp_points_face_"<<std::endl;
//     	// std::cout<<reference_qp_points_face_<<std::endl;
//         shape_face_.init(coeff_,reference_qp_points_face_);
//         // std::cout<<"shape_face_.eval()"<<std::endl;
//         // std::cout<<shape_face_.eval()<<std::endl;
//         // Vector<Real,NDofsFace> out; 

//         const auto& normal=signed_normal.normals()[F_FE.elem_id()][F_FE.side_id()];
        
//         // std::cout<<"qui1"<<std::endl;
//         Integer id=F_FE.side_id();
//         // std::cout<<"qui2"<<std::endl;
//         trace_map_.init(F_FE);
//         // std::cout<<"qui3"<<std::endl;
//         shape_trace_face_.init_map(trace_map_);
//         // std::cout<<"qui4"<<std::endl;
//         shape_trace_face_.init(coeff_,id);
//         // std::cout<<"F_FE.side_id()"<<std::endl;
//         // std::cout<<F_FE.side_id()<<std::endl;

//         // std::cout<<"shape_trace_face_.eval()"<<std::endl;
//         // std::cout<<shape_trace_face_.eval()<<std::endl;
        
//         for(Integer i=0;i<dofs_per_entity[0]; i++)
//         {
//          compute_face_dofs(mat,cont,QRuleFace::qp_weights,shape_face_.eval() ,shape_trace_face_.eval()[i], normal, F_FE);
//          cont++;
//         }
//     	// std::cout<<"mat"<<std::endl;
//     	// std::cout<<mat<<std::endl;

 


//      // qp_points_face_
//      // shape_face_.init_map(C_map);
//     }


//     // momentum dofs
//     template<Integer EntityDim, typename Mat>
//     std::enable_if_t< (EntityDim==entity[1]), void> 
//     compute(Mat& mat, Integer& cont, const FiniteElem<Elem>& C_FE, FiniteElem<Elem>& F_FE)
//     {
//      // qp_points_momentum_
//      auto mesh_ptr=F_FE.mesh_ptr();
//      const auto& signed_normal=mesh_ptr->signed_normal();
//      C_map_.init(C_FE);
//      F_Pn_map_.init(F_FE);
//      shape_momentum_.init_map(C_map_);
//      shape_lagrange_momentum_.init_map(F_Pn_map_);
//      SingleShapeFunctionCoefficientsCollection<Elem,FEFamily,Order>::apply(signed_normal.sign(C_FE.elem_id()),coeff_);
//      F_FE.transform_point(qp_points_momentum_,QRuleMomentum::qp_points);
//      C_FE.points_to_reference_points(reference_qp_points_momentum_,qp_points_momentum_);

//      // std::cout<<"QRuleMomentum::qp_points"<<std::endl;
//      // std::cout<<QRuleMomentum::qp_points<<std::endl;
//      // std::cout<<"qp_points_momentum_"<<std::endl;
//      // std::cout<<qp_points_momentum_<<std::endl;
//      // std::cout<<"reference_qp_points_momentum_"<<std::endl;
//      // std::cout<<reference_qp_points_momentum_<<std::endl;

//      shape_momentum_.init(reference_qp_points_momentum_);
//      shape_lagrange_momentum_.init(reference_qp_points_momentum_);
     
//      int cont2=0;
//      // Vector<Real,NDofsMomentum > out; 
//      for(Integer i=NDofsFace;i<NDofs;i++)
//      {
//      compute_momentum_dofs(mat,cont,QRuleFace::qp_weights,shape_momentum_.eval(),cont2,shape_lagrange_momentum_.eval()[cont2], F_FE);
//      cont++;
//      cont2++;

//      }


//      // std::cout<<"QRuleMomentum::qp_points"<<std::endl;
//      // std::cout<<QRuleMomentum::qp_points<<std::endl;

//      // std::cout<<"qp_points_momentum_"<<std::endl;
//      // std::cout<<qp_points_momentum_<<std::endl;



//      // std::cout<<"shape_momentum_.eval()"<<std::endl;
//      // std::cout<<shape_momentum_.eval()<<std::endl;
//      // std::cout<<"shape_lagrange_momentum_.eval()"<<std::endl;
//      // std::cout<<shape_lagrange_momentum_.eval()<<std::endl;

//      // std::cout<<"mat"<<std::endl;
//      // std::cout<<mat<<std::endl;
//     }

//     // template<Integer H,typename Mat>
//     // std::enable_if_t<H==NDofs,void> 
//     // compute_aux(Mat& mat, Integer& cont, const FiniteElem<Elem>& C_FE, FiniteElem<Elem>& F_FE)
//     // {} 


//     // template<Integer H,typename Mat>
//     // std::enable_if_t<H<NDofs,void> 
//     // compute_aux(Mat& mat, Integer& cont, const FiniteElem<Elem>& C_FE, FiniteElem<Elem>& F_FE)
//     // {
//     // 	compute<entity[0]>(mat,cont,C_FE,F_FE);
//     // 	compute<entity[1]>(mat,cont,C_FE,F_FE);
//     // 	cont++;
//     // 	compute_aux<H+1>(mat,cont,C_FE,F_FE);
//     // } 


//     template<typename Mat>
//     void compute(Mat& mat, const FiniteElem<Elem>& C_FE, FiniteElem<Elem>& F_FE)
//     {
//     	Integer cont=0;
//     	// compute_aux<Cont>(mat,cont,C_FE,F_FE);
//     	// compute_face_aux<0>(mat,cont,C_FE,F_FE);
//     	// compute_momentum_aux<0>(mat,cont,C_FE,F_FE);
//         // std::cout<<"compute first"<<std::endl;

//         compute<entity[1]>(mat,cont,C_FE,F_FE);

//     	for(Integer i=0;i<ManifoldDim+1;i++)
//     		{F_FE.init_boundary(i,true);
//     	     compute<entity[0]>(mat,cont,C_FE,F_FE);}

    	


//      // std::cout<<"mat"<<std::endl;
//      // std::cout<<mat<<std::endl;

//     	// compute<entity[1]>(mat,C_FE,F_FE);
//     }  


// private:
// 	ShapeFace shape_face_;
// 	ShapeTraceFace shape_trace_face_;
// 	ShapeMomentum shape_momentum_;
// 	ShapeLagrangianMomentum shape_lagrange_momentum_;
// 	QPpointsFace qp_points_face_;
// 	QPpointsFace reference_qp_points_face_;
// 	QPpointsMomentum qp_points_momentum_;
// 	QPpointsMomentum reference_qp_points_momentum_;
// 	Map C_map_;
// 	TraceMap trace_map_;
// 	Map F_map_;
// 	F_Lagrangian_Map F_Pn_map_;
// 	// Coefficients coeffs_;
// 	Array<Real,NDofs> coeff_;
// };



template <Integer Dim,Integer ManifoldDim, Integer FEFamily, Integer Order,Integer Continuity, Integer NComponents_,Integer QR>
class Dof<ElementFunctionSpace<Simplex<Dim,ManifoldDim>,FEFamily,Order,Continuity,NComponents_>,QR>
{
 public:
 	// static constexpr Integer FEFamily=LagrangeFE;
 	using Elem=Simplex<Dim,ManifoldDim>;



 	// HA SENSO PER NCOMPONENTS == 1 RIFLETTI PER PIU COMPONENTI

 	static constexpr Integer NComponents=1; 

 	using BaseFunctionSpace=BaseFunctionSpace<FEFamily,Order,Continuity,NComponents>;
 	using ElementFunctionSpace=ElementFunctionSpace<Elem,FEFamily,Order,Continuity,NComponents>;
 	static constexpr auto entity=ElementFunctionSpace::entity;
    static constexpr Integer QRuleOrder=QuadratureOrder<IdentityOperator, BaseFunctionSpace>::value;
    using QRule=typename QuadratureRule<QR>:: template rule<Elem,1>;
    using Shape=ShapeFunction<Elem,BaseFunctionSpace,IdentityOperator,QRule>;
    using ShapeValue=typename Shape::type;
    using Map=MapFromReference<IdentityOperator,Elem,FEFamily> ;
    static constexpr auto Dofs=LagrangeDofs<ElementFunctionSpace>::value();
    static constexpr auto NDofs=LagrangeDofs<ElementFunctionSpace>::NDofs;
    
    template<typename MeshT>
    Dof(const MeshT& mesh):
    C_FE_(mesh),
    F_FE_(mesh)
    {}

    template<Integer M,Integer N=0>
    inline std::enable_if_t<N==NDofs,void>
    loop_coarse_dof(const Elem& C_elem, const Elem& F_elem)
    {}


    template<Integer M,Integer N=0>
    inline std::enable_if_t<N<NDofs,void>
    loop_coarse_dof(const Elem& C_elem, const Elem& F_elem)
    {
     // HA SENSO PER NCOMPONENTS == 1 RIFLETTI PER PIU COMPONENTI
     // TODO FIXME
     // auto& value=func_values_[N];

     if(FEFamily==LagrangeFE)
     	local_mat_(M,N)=func_values_[N][0](0,0);
     // DofAux<M,N,ElementFunctionSpace,QR>::compute(local_mat_,func_values_);



     // std::cout<<"loop_coarse_dof =="<<N<<std::endl;
     // std::cout<<value<<std::endl;
     // local_mat_(M,N)=func_values_[N][0](0,0);
     // shape_.init();
     // loop on coarse shape functions
     loop_coarse_dof<M,N+1>(C_elem,F_elem);
        
    }





    
    void compute_point_row()
    {
    const auto& J_inv=C_FE_.inv_jac();
    const auto& v0=C_FE_.v0();
    const auto& rows=J_inv.rows();
    const auto& cols=J_inv.cols();

    for(Integer i=0;i<rows;i++)
    	{
    		transformed_point_row_(0,i)=J_inv(i,0)*(transformed_point_col_[0]-v0[0]);
    	}

    for(Integer i=0;i<rows;i++)
    	{
	    	for(Integer j=1;j<cols;j++)
	    	{
	    		transformed_point_row_(0,i)+=J_inv(i,j)*(transformed_point_col_[j]-v0[j]);
	    	}
       	}

     }



    template<Integer N=0>
    inline std::enable_if_t<N==NDofs,void>
    loop_fine_dof(const Elem& C_elem, const Elem& F_elem)
    {}

    template<Integer N=0>
    inline std::enable_if_t<N<NDofs,void>
    loop_fine_dof(const Elem& C_elem, const Elem& F_elem)
    {
     std::cout<<"loop_fine_dof =="<<N<<"    "<<Dofs[N]<<std::endl;
     F_FE_.init(F_elem.id);
     F_map_.init(F_FE_);
     F_FE_.transform_point(transformed_point_col_,Dofs[N]);
     // transformed_point_row_=C_FE_.inv_jac()*transformed_point_col_;
     compute_point_row();

     shape_.init(transformed_point_row_,C_FE_);
     // std::cout<<"Dofs[N] =="<<std::endl;
     // std::cout<<Dofs[N]<<std::endl;
     // std::cout<<"transformed_point_col_ =="<<std::endl;
     // std::cout<<transformed_point_col_ <<std::endl;
     // std::cout<<"transformed_point_row_ =="<<std::endl;
     // std::cout<<transformed_point_row_ <<std::endl;
     func_values_=shape_.eval();
     // std::cout<<"func_values_ =="<<std::endl;
     // std::cout<<func_values_<<std::endl;
     loop_coarse_dof<N>(C_elem,F_elem);
     // loop on coarse shape functions
     loop_fine_dof<N+1>(C_elem,F_elem);
        
    }


    inline void compute(const Elem& C_elem, const Elem& F_elem)
    {

     C_FE_.init(C_elem.id);
     C_FE_.init_inv_jac();
     C_map_.init(C_FE_);
     shape_.init_map(C_map_);




     // loop on coarse shape functions
     // std::cout<<"compute"<<std::endl;
     // std::cout<<"map_"<<C_map_()<<std::endl;
     // std::cout<<"Dofs reference =="<<std::endl;
     // std::cout<<Dofs<<std::endl;
     loop_fine_dof(C_elem,F_elem);
     std::cout<<"local_mat_"<<std::endl;
     std::cout<<local_mat_<<std::endl;
    }

 private:
 	Shape shape_;
 	Map C_map_;
 	Map F_map_;
    FiniteElem<Elem> C_FE_;
    FiniteElem<Elem> F_FE_;
    ShapeValue func_values_;
    Matrix<Real,1,Dim> transformed_point_row_;
    Matrix<Real,Dim,1> transformed_point_col_;
    Matrix<Real,NDofs,NDofs> local_mat_;
};

























template<typename FunctionSpace>
class FullFunctionSpaceInterpolation
{
public:
	using Elem=typename FunctionSpace::Elem;
	using TupleOfSpaces=typename FunctionSpace::FunctionSpace::TupleOfSpaces;
	using ElemDofMap=typename FunctionSpace::DofsDM::ElemDofMap;
    static constexpr Integer  NLocalDofs=FunctionSpace::DofsDM::NLocalDofs;
	FullFunctionSpaceInterpolation(std::shared_ptr<FunctionSpace> W_ptr):
	spaces_ptr_(W_ptr),
	C_FE_(W_ptr->mesh()),
	F_FE_(W_ptr->mesh())
	{}

    inline void add_functionspace(std::shared_ptr<FunctionSpace> W_ptr){spaces_ptr_=W_ptr;}


   template<Integer N=0, Integer Nmax=TupleTypeSize<TupleOfSpaces>::value>
   inline std::enable_if_t<N==Nmax,void> 
   space_loop(SparseMatrix<Real>& A, const Integer C_el, const Integer C_level, const Integer F_el,const Integer F_level)
   {}


   template<Integer N=0, Integer Nmax=TupleTypeSize<TupleOfSpaces>::value>
   inline std::enable_if_t<N<Nmax,void> 
   space_loop(SparseMatrix<Real>& A, const Integer C_el, const Integer C_level, const Integer F_el,const Integer F_level)
   {
     using Space=GetType<TupleOfSpaces,N>;
     constexpr Integer n_dofs=FunctionSpaceDofsPerElem<Space>::value;
     constexpr Integer NComponents=Space::NComponents;
     constexpr Integer n_dofs_one_component=n_dofs/NComponents;
     Matrix<Real,n_dofs_one_component,n_dofs_one_component> mat;
   	 DofAux<Space> dof;
   	 

   	 C_FE_.init(C_el,C_level);
   	 F_FE_.init(F_el,F_level);

   	 auto& dm_row=tuple_get<N>(elemdofmap_row_);
   	 auto& dm_col=tuple_get<N>(elemdofmap_col_);


   	 auto& dofsdofmap=spaces_ptr_->dofsdofmap();
   	 
   	 // std::cout<< " Space ==  " << N << std::endl;
   	 // std::cout<< " F_FE_.elem_id() ==  " << F_FE_.elem_id() << std::endl;
   	 // std::cout<< " C_FE_.elem_id() ==  " << C_FE_.elem_id() << std::endl;
   	 // std::cout<< " F_FE_.level() ==  " << F_FE_.level() << std::endl;
   	 // std::cout<< " C_FE_.level() ==  " << C_FE_.level() << std::endl;
   	 // std::cout<< " n_dofs ==  " << n_dofs << std::endl;
     
   	 dofsdofmap.template dofmap_get<N>(dm_row,F_FE_.elem_id(),F_FE_.level());
   	 dofsdofmap.template dofmap_get<N>(dm_col,C_FE_.elem_id(),C_FE_.level());


     F_FE_.init(F_FE_.elem_id(),F_FE_.level());
     // std::cout<< "pre dof compute " <<std::endl;
     dof.compute(mat,C_FE_,F_FE_);
     // std::cout<< "after dof compute " <<std::endl;
     // std::cout<< "mat " <<std::endl;
     // std::cout<< mat  <<std::endl;
     // std::cout<< "n_dofs " <<std::endl;
     // std::cout<< n_dofs  <<std::endl;
     // std::cout<< "n_dofs_one_component " <<std::endl;
     // std::cout<< n_dofs_one_component  <<std::endl;
     // std::cout<< "NComponents " <<std::endl;
     // std::cout<< NComponents  <<std::endl;
     // std::cout<< "dm_row" <<std::endl;
     // std::cout<< dm_row <<std::endl;
     // std::cout<< "dm_col" <<std::endl;
     // std::cout<< dm_col <<std::endl;
         
     for(Integer i=0;i<n_dofs_one_component;i++)
     	for(Integer j=0;j<n_dofs_one_component;j++)
     		for(Integer k=0;k<NComponents;k++)
     		{
              // std::cout<<"("<<dm_row[i*NComponents+k]<<", "<<dm_col[j*NComponents+k]<<", "<<mat(i,j)<<")   ";
              A.equal(mat(i,j),dm_row[i*NComponents+k],dm_col[j*NComponents+k]);
              // A.print_row(dm_row[i*NComponents+k]);
     		}
     		
     // A.print_val();
     space_loop<N+1,Nmax>(A,C_el,C_level,F_el,F_level);
   }

    
    
    void find_children(SparseMatrix<Real>& A,const Integer C_level,const Integer F_level, const Integer el,const Integer old_el)
    {
        auto &mesh=spaces_ptr_->mesh();
        auto &bisection=spaces_ptr_->bisection();
        auto &tracker=bisection.tracker();
        // std::cout<<"F_level=="<<F_level<<std::endl;
        // std::cout<<"find_children el=="<<el<<std::endl;
    	auto& fine_elem=mesh.elem(el);
      	auto& children=fine_elem.children;
      	// for(std::size_t i=0;i<children.size();i++)
      	// 	std::cout<<children[i]<<" ";
      	// std::cout<<std::endl;

        if(elem_belongs_to_level(mesh,el,F_level,tracker))
            {
            F_FE_.init(el,F_level);
      		space_loop(A,old_el,C_level,el,F_level);
            }
        else
        	{
	      	for(std::size_t i=0;i<children.size();i++)
	      	{
	            
	      		if(elem_belongs_to_level(mesh,children[i],F_level,tracker))
	      		{
	               
	      		   space_loop(A,old_el,C_level,children[i],F_level);
	      		}
	      		else
	      		{
	      			find_children(A,C_level,F_level,children[i],old_el);
	      		}


	      	}
        	}

    }


	void init(SparseMatrix<Real> &A,const Integer C_level, const Integer F_level)
	{
      

      // SparseMatrix<Real> A;
      auto &mesh=spaces_ptr_->mesh();
      auto &bisection=spaces_ptr_->bisection();
      auto &tracker=bisection.tracker();
      auto& dofsdofmap=spaces_ptr_->dofsdofmap();
      auto& level_cumultive_n_dofs=dofsdofmap.level_cumultive_n_dofs();
      auto& level_cumulative_dofs_array=dofsdofmap.level_cumulative_dofs_array();
      auto& n2e=spaces_ptr_->node2elem();
      // std::cout<<" qui1 "<<std::endl;
      Integer n_dofs_rows=level_cumultive_n_dofs[F_level];
      // std::cout<<" qui2 "<<std::endl;
      Integer n_dofs_cols=level_cumultive_n_dofs[C_level];
      // std::cout<<" qui3 ,n_dofs_rows= "<<n_dofs_rows<<std::endl;
      Integer max_cols=n2e.max_n_nodes();
      // std::cout<<" qui4 , n_dofs_cols ="<<n_dofs_cols<<std::endl;
      max_cols=min(NLocalDofs*max_cols,n_dofs_rows);
      // std::cout<<" qui5 ,max_cols ="<<max_cols<<std::endl;
      level_cumultive_n_dofs[C_level];
      
      // std::cout<<" qui6 "<<std::endl;
      A.init(n_dofs_rows,n_dofs_cols,max_cols);
      // std::cout<<"n_dofs_rows = "<<n_dofs_rows<<std::endl;
      // std::cout<<"max_cols = "<<max_cols<<std::endl;


      // std::cout<<"level_cumulative_dofs_array"<<std::endl;
      // for(Integer i=0;i<level_cumulative_dofs_array.size();i++)
      // 	{for(Integer j=0;j<level_cumulative_dofs_array[i].size();j++)
      //    std::cout<<level_cumulative_dofs_array[i][j]<<" ";
      //    std::cout<<std::endl;}
      for(std::size_t C_el=0;C_el<mesh.n_elements();C_el++)
      {
      	// std::cout<<"C_el=="<<C_el<<std::endl;
      	// if the elem does not belong to the C_level, continue

      	// if(C_el==132)
      	// {
      
       //  std::cout<<"START BUG C_el=="<<C_el<<std::endl;
       //  std::cout<<"C_level=="<<C_level<<std::endl;
       //  std::cout<<"track.get_iterate(i)=="<<tracker.get_iterate(C_el)<<std::endl;
       
      	if(!elem_belongs_to_level(mesh,C_el,C_level,tracker))continue;
        // std::cout<<"eeeeeevviva C_el=="<<C_el<<std::endl;
        C_FE_.init(C_el,C_level);
        C_FE_.init_inv_jac();
        // auto& elem=mesh.elem(el);
        find_children(A,C_level,F_level,C_el,C_el);
     	// std::cout<<"END BUG C_el=="<<C_el<<std::endl;
      	// }

      	// if(!elem_belongs_to_level(mesh,C_el,C_level,tracker))continue;
       //  std::cout<<"C_el=="<<C_el<<std::endl;
       //  C_FE_.init(C_el,C_level);
       //  C_FE_.init_inv_jac();
       //  find_children(A,C_level,F_level,C_el,C_el);




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
 FiniteElem<Elem> C_FE_;
 FiniteElem<Elem> F_FE_;
 ElemDofMap elemdofmap_row_;
 ElemDofMap elemdofmap_col_;

};

template<typename FunctionSpace>
class FullFunctionSpaceLevelsInterpolation
{
public:
	using FunctionSpaceInterpolation=FullFunctionSpaceInterpolation<FunctionSpace>;
    FullFunctionSpaceLevelsInterpolation(std::shared_ptr<FunctionSpace> W_ptr):
	spaces_ptr_(W_ptr)
	{}

	// levels contain


    inline auto& matrices(){return mat_vec_;}

    inline auto& matrix(const Integer i){return mat_vec_[i];}

    inline auto& matrix(const Integer i)const{return mat_vec_[i];}

    inline auto& matrix(const Integer coarse,const Integer fine)
    {   for(Integer i=0;i<interp_levels_vec_.size();i++)
    	{
    	if(interp_levels_vec_[i][0]==coarse && interp_levels_vec_[i][1]==fine)
    		return mat_vec_[i];
    	}
    	std::cout<<"ERROR: not existing coarse-fine levels pair in FullFunctionSpaceLevelsInterpolation"<<std::endl;
    }

	void init(const std::vector<Integer> levels)
	{
        Integer n_levels_=levels.size();
        mat_vec_.resize(n_levels_-1);
        interp_levels_vec_.resize(n_levels_-1);
        levels_.resize(n_levels_);

        for(Integer i=0;i<n_levels_;i++)
        	levels_[i]=levels[i];

		for(Integer i=0;i<n_levels_-1;i++)
		{
			interp_vec_.push_back(FunctionSpaceInterpolation(spaces_ptr_));
			std::cout << " FullFunctionSpaceLevelsInterpolation " <<i<< std::endl;
			interp_levels_vec_[i][0]=levels[i+1];
			interp_levels_vec_[i][1]=levels[i+1];
			interp_vec_[i].init(mat_vec_[i],levels[i],levels[i+1]);

		}
	}

    
    // auto build_single_interp(SparseMatrix<Real>& A_C,const SparseMatrix<Real>& A_F, const Integer level)
    // {
     
    //  A_C


    // }








private:
	Integer n_levels_;
	std::vector<Integer> levels_;
	std::shared_ptr<FunctionSpace> spaces_ptr_;
	std::vector<FunctionSpaceInterpolation> interp_vec_;
	std::vector<SparseMatrix<Real>> mat_vec_;
	std::vector<Array<Integer,2>> interp_levels_vec_;
};




   void plus_equal(std::vector<Real>& x, const std::vector<Real>& b)
   {
   	for(Integer i=0;i<x.size();i++)
   		x[i]+=b[i];
   }


  Real l2_norm(const std::vector<Real>& vec)
  {
  	Real tmp=0.0;

  	for(Integer i=0;i<vec.size();i++)
  		tmp+=(vec[i]*vec[i]);

    return sqrt(tmp);

  }


   template<typename FunctionSpace>
   inline void vcycle(std::vector<Real>& x,
   	                     const std::vector<SparseMatrix<Real>>& A_levels,
   	                     const std::vector<Real>& b,
   	                     const FullFunctionSpaceLevelsInterpolation<FunctionSpace>& interpolation,
   	                     const Integer pre_smoothing,
   	                     const Integer post_smoothing,
   	                     const Integer level)
   {
    std::vector<Real> F_rhs;
    std::vector<Real> C_rhs;
    std::vector<Real> C_correction;
    std::vector<Real> F_correction;
    // std::cout<<"vcycle level = "<<level <<std::endl;

   	if(level==0)
   	{
   		// std::cout<<"vcycle coarse level = "<<level <<std::endl;
     gauss_seidel(x,A_levels[0],b,10000);
   	}
   	else
   	{

   	//  std::cout<<"vcycle pre smoothing " <<std::endl;
    // for(Integer i=0;i<x.size();i++)
   	//  	std::cout<<x[i]<<std::endl;
   	 gauss_seidel(x,A_levels[level],b,pre_smoothing);
   	//  std::cout<<"vcycle post pre smoothing " <<std::endl;
    // for(Integer i=0;i<x.size();i++)
   	//  	std::cout<<x[i]<<std::endl;
   	 // std::cout<<"vcycle pre rhs " <<std::endl;
   	 A_levels[level].multiply_and_add(F_rhs,-1.0,x,b);
   	 interpolation.matrix(level-1).transpose_and_multiply(C_rhs,F_rhs);
   	 // std::cout<<"vcycle go from"<<level<<" to "<<level -1<<std::endl;
   	 // x_interp_vec[i]=interpolation.matrix(i-1).multiply(x_interp_vec[i-1]);


     vcycle(C_correction,A_levels,C_rhs,interpolation,pre_smoothing,post_smoothing,level-1);
     // std::cout<<"vcycle interpolation "<<std::endl;

     interpolation.matrix(level-1).multiply(F_correction,C_correction);
   	 // std::cout<<"C_correction " <<std::endl;
     // for(Integer i=0;i<C_correction.size();i++)
   	 // 	std::cout<<C_correction[i]<<std::endl;
   	 // std::cout<<"F_correction " <<std::endl;
     // for(Integer i=0;i<F_correction.size();i++)
   	 // 	std::cout<<F_correction[i]<<std::endl;
   	 // std::cout<<"x pre plus_equal  " <<std::endl;
     // for(Integer i=0;i<x.size();i++)
   	 // 	std::cout<<x[i]<<std::endl;
     plus_equal(x,F_correction);
   	 // std::cout<<"vcycle post plus_equal " <<std::endl;
     // for(Integer i=0;i<x.size();i++)
   	 // 	std::cout<<x[i]<<std::endl;

     // std::cout<<"vcycle pre post smoothing "<<std::endl;
     gauss_seidel(x,A_levels[level],b,post_smoothing);
   	 // std::cout<<"vcycle post post smoothing " <<std::endl;
     // for(Integer i=0;i<x.size();i++)
   	 // 	std::cout<<x[i]<<std::endl;
   	}

   }

  template<typename FunctionSpace,typename VecVecVec>
   inline void patch_vcycle(std::vector<Real>& x,
   	                     const std::vector<SparseMatrix<Real>>& A_levels,
   	                     const std::vector<Real>& b,
   	                     const FullFunctionSpaceLevelsInterpolation<FunctionSpace>& interpolation,
   	                     const VecVecVec& nodes2entity,
   	                     const Integer pre_smoothing,
   	                     const Integer post_smoothing,
   	                     const Integer level)
   {
    std::vector<Real> F_rhs;
    std::vector<Real> C_rhs;
    std::vector<Real> C_correction;
    std::vector<Real> F_correction;
    // std::cout<<"vcycle level = "<<level <<std::endl;

   	if(level==0)
   	{
   		// std::cout<<"vcycle coarse level = "<<level <<std::endl;
     // gauss_seidel(x,A_levels[0],b,10000);
   		patch_gauss_seidel(x,A_levels[0],b,nodes2entity[0],100);
   	}
   	else
   	{

   	//  std::cout<<"vcycle pre smoothing " <<std::endl;
    // for(Integer i=0;i<x.size();i++)
   	//  	std::cout<<x[i]<<std::endl;
   	 patch_gauss_seidel(x,A_levels[level],b,nodes2entity[level],pre_smoothing);
   	//  std::cout<<"vcycle post pre smoothing " <<std::endl;
    // for(Integer i=0;i<x.size();i++)
   	//  	std::cout<<x[i]<<std::endl;
   	 // std::cout<<"vcycle pre rhs " <<std::endl;
   	 A_levels[level].multiply_and_add(F_rhs,-1.0,x,b);
   	 interpolation.matrix(level-1).transpose_and_multiply(C_rhs,F_rhs);
   	 // std::cout<<"vcycle go from"<<level<<" to "<<level -1<<std::endl;
   	 // x_interp_vec[i]=interpolation.matrix(i-1).multiply(x_interp_vec[i-1]);

   	 C_correction.resize(C_rhs.size(),0.0);


     patch_vcycle(C_correction,A_levels,C_rhs,interpolation,nodes2entity,pre_smoothing,post_smoothing,level-1);
     interpolation.matrix(level-1).multiply(F_correction,C_correction);
   	 plus_equal(x,F_correction);
    

     // std::cout<<"vcycle interpolation "<<std::endl;





   	 // std::cout<<"C_correction " <<std::endl;
     // for(Integer i=0;i<C_correction.size();i++)
   	 // 	std::cout<<C_correction[i]<<std::endl;
   	 // std::cout<<"F_correction " <<std::endl;
     // for(Integer i=0;i<F_correction.size();i++)
   	 // 	std::cout<<F_correction[i]<<std::endl;
   	 // std::cout<<"x pre plus_equal  " <<std::endl;
     // for(Integer i=0;i<x.size();i++)
   	 // 	std::cout<<x[i]<<std::endl;
     
   	 // std::cout<<"vcycle post plus_equal " <<std::endl;
     // for(Integer i=0;i<x.size();i++)
   	 // 	std::cout<<x[i]<<std::endl;

     // std::cout<<"vcycle pre post smoothing "<<std::endl;
     patch_gauss_seidel(x,A_levels[level],b,nodes2entity[level],post_smoothing);
   	 // std::cout<<"vcycle post post smoothing " <<std::endl;
     // for(Integer i=0;i<x.size();i++)
   	 // 	std::cout<<x[i]<<std::endl;
   	}

   }


   template<typename FunctionSpace>
   inline void multigrid(std::vector<Real>& x,
   	                     const std::vector<SparseMatrix<Real>>& A_levels,
   	                     const std::vector<Real>& b,
   	                     const FullFunctionSpaceLevelsInterpolation<FunctionSpace>& interpolation,
   	                     const Integer pre_smoothing,
   	                     const Integer post_smoothing,
   	                     const Integer level,
   	                     const Integer max_iter,
   	                     const Real toll)
   {
   	std::vector<Real> rhs;
   	Real norm_rhs;
   	// std::cout<<"x size"<<x.size()<<std::endl;
   	// std::cout<<"b size"<<b.size()<<std::endl;
   	// std::cout<<"A rows"<<A_levels[level].max_rows()<<std::endl;
   	// std::cout<<"A cols"<<A_levels[level].max_cols()<<std::endl;

   	// std::cout<<"print b in multigrid" <<std::endl;
    // for(Integer i=0;i<x.size();i++)
   	//  	std::cout<<b[i]<<std::endl;
     
   	for(Integer i=0;i<max_iter;i++)
   	{
   	
     A_levels[level].multiply_and_add(rhs,-1.0,x,b);
     norm_rhs=l2_norm(rhs);
     std::cout<<"multigrid iteration == "<< i<<"  pre norm_rhs="<<norm_rhs<<std::endl;
     // std::cout<<"pre_smoothing, post_smoothing"<< pre_smoothing<<" , "<< post_smoothing<<std::endl;
   	 vcycle(x,A_levels,b,interpolation,pre_smoothing,post_smoothing,level);
    // std::cout<<"print x in multigrid" <<std::endl;
    // for(Integer i=0;i<x.size();i++)
   	//  	std::cout<<x[i]<<std::endl;

   	 // A_levels[level].print_val();

     A_levels[level].multiply_and_add(rhs,-1.0,x,b);
     norm_rhs=l2_norm(rhs);
     std::cout<<"multigrid iteration == "<< i<<"  after norm_rhs="<<norm_rhs<<std::endl;
     if(norm_rhs<toll)
     	break;
   	}

   }








   template<typename FunctionSpace,typename VecVecVec>
   inline void patch_multigrid(std::vector<Real>& x,
   	                     const std::vector<SparseMatrix<Real>>& A_levels,
   	                     const std::vector<Real>& b,
   	                     const FullFunctionSpaceLevelsInterpolation<FunctionSpace>& interpolation,
   	                     const VecVecVec& entity2dofs,
   	                     const Integer pre_smoothing,
   	                     const Integer post_smoothing,
   	                     const Integer level,
   	                     const Integer max_iter,
   	                     const Real toll)
   {
   	std::vector<Real> rhs;
   	Real norm_rhs;
   	// std::cout<<"x size"<<x.size()<<std::endl;
   	// std::cout<<"b size"<<b.size()<<std::endl;
   	// std::cout<<"A rows"<<A_levels[level].max_rows()<<std::endl;
   	// std::cout<<"A cols"<<A_levels[level].max_cols()<<std::endl;

   	// std::cout<<"print b in multigrid" <<std::endl;
    // for(Integer i=0;i<x.size();i++)
   	//  	std::cout<<b[i]<<std::endl;
     std::vector<Real> norm_rhs_vec;
   	for(Integer i=0;i<max_iter;i++)
   	{
   	
     A_levels[level].multiply_and_add(rhs,-1.0,x,b);
     norm_rhs=l2_norm(rhs);
     std::cout<<"PATCH multigrid iteration == "<< i<<"  pre norm_rhs="<<norm_rhs<<std::endl;

   	 patch_vcycle(x,A_levels,b,interpolation,entity2dofs,pre_smoothing,post_smoothing,level);


     A_levels[level].multiply_and_add(rhs,-1.0,x,b);
     norm_rhs=l2_norm(rhs);
     norm_rhs_vec.push_back(norm_rhs);
     std::cout<<"PATCH multigrid iteration == "<< i<<"  after norm_rhs="<<norm_rhs<<std::endl;
     if(norm_rhs<toll)
     	break;
   	}
   	for(Integer i=0;i<norm_rhs_vec.size();i++)
   		std::cout<<norm_rhs_vec[i]<<std::endl;

   }



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





	class ExactLinearElasticity2D
	{
	public: 
	    // using Point=Matrix<Real,3,1>;
		using type=Matrix<Real,2,1>;
	    template<typename Point>
		static auto eval(const Point& p)
		{
			// auto M_PI_Squared=M_PI*M_PI;
			const auto& x=p[0];
			const auto& y=p[1];

			type func{(1.0)*( 4.0 * M_PI_Squared * sin(M_PI * x) * sin(M_PI * y) - 2.0 * M_PI_Squared *cos( M_PI * x) * cos( M_PI * y) ),
                      (1.0)*( 4.0 * M_PI_Squared * sin(M_PI * x) * sin(M_PI * y) - 2.0 * M_PI_Squared *cos( M_PI * x) * cos( M_PI * y) )
			          };
			return func; 
		}
	};

	class ExactLinearElasticity3D
	{
	public: 
	    // using Point=Matrix<Real,3,1>;
		using type=Matrix<Real,3,1>;


	    template<typename Point,typename FiniteElem>
		static type eval(const Point& p,FiniteElem& FE)
		{
			const auto& x=p[0];
			const auto& y=p[1];
			const auto& z=p[2];

			Real A=5.0*M_PI_Squared*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
			Real B=2.0*M_PI_Squared*cos(M_PI*y)*cos(M_PI*z)*sin(M_PI*x);
			Real C=2.0*M_PI_Squared*cos(M_PI*x)*cos(M_PI*z)*sin(M_PI*y);

			type func{(A - B - C ),(A - B - C ),(A - B - C )};
				
			return func; 
		}
	};


	template<Integer ManifoldDim>
	class ExactLinearElasticity
	{
	public: 
		using type=Matrix<Real,ManifoldDim,1>;

	    template<Integer N,typename Point>
		static std::enable_if_t<N==2,type> eval_aux(const Point& p)
		{
				return ExactLinearElasticity2D::eval(p);  
		}

	    template<Integer N,typename Point>
		static std::enable_if_t<N==3,type> eval_aux(const Point& p)
		{
				return ExactLinearElasticity3D::eval(p);  
		}

	    template<typename Point,typename FiniteElem>
		static type eval(const Point& p,FiniteElem& FE)
		{
			return eval_aux<ManifoldDim>(p);
		}

	};




template<Integer Dim>
class NormalFunction;

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
	void LSFEM_ContactLinearElasticity(const Integer n, const Integer level, const Integer n_levels)
	{
	  // constexpr Integer ManifoldDim=2;
		constexpr Integer Dim=ManifoldDim;
		using MeshT=Mesh<Dim, ManifoldDim>;
		using Elem=typename MeshT::Elem;
		MeshT mesh;
		if(ManifoldDim==2)
		{
			std::cout<<"2D LSFEM_ContactLinearElasticity";
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
			std::cout<<"3D LSFEM_ContactLinearElasticity";
			generate_cube(mesh,n,n,n);
			mesh.build_dual_graph();
			mark_boundary(mesh);

	      // read_mesh("../data/tetrahedron.MFEM", mesh);
	      // assign_tags(mesh);
	      // Bisection<MeshT> bisection(mesh);
	      // bisection.uniform_refine(0);
	      // mesh.update_dual_graph();
		}
        // Integer n_levels=12;
	    clock_t begin=clock();
		clock_t end = clock();
		double elapsed_secs;

		// std::cout<<"mat==="<<mat3<<std::endl;

		

        std::cout<<"n2em init="<<std::endl;


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

	    Node2ElemMap<MeshT> n2em(mesh,bisection);
	    n2em.init();	      
	    // n2em.add_bisection(bisection);
	    // n2em.init();
	    std::cout<<"--------------:mesh n elements=="<<mesh.n_elements()<<std::endl;
	    std::cout<<"--------------:n2em n max_n_nodes=="<<n2em.max_n_nodes()<<std::endl;






		using AuxRT_n= FunctionSpace< MeshT, RT<Order1,ManifoldDim>>;
		using AuxP_n= FunctionSpace< MeshT, Lagrange<Order2,ManifoldDim>>;
		using AuxP_0= FunctionSpace< MeshT, Lagrange<0,ManifoldDim>>;
		using LSFEM= FunctionSpace< MeshT, RT<Order1,ManifoldDim>,Lagrange<Order2,ManifoldDim>>;

		AuxRT_n rtn(mesh,bisection,n2em);//csmc);
		AuxP_n pn(mesh,bisection,n2em);//csmc);
		AuxP_0 p0(mesh,bisection,n2em);//csmc);
		std::cout<<"FIRST PRE UPDATE="<<std::endl;
		LSFEM lsfem(mesh,bisection,n2em);//csmc);
		std::cout<<"FIRST POST UPDATE="<<std::endl;
	    std::cout<<"--------------:rtn n max_n_nodes=="<<rtn.node2elem().max_n_nodes()<<std::endl;
	    std::cout<<"--------------:pn n max_n_nodes=="<<pn.node2elem().max_n_nodes()<<std::endl;

	    auto Wtrial=MixedFunctionSpace(rtn,pn);
		std::cout<<"FIRST POST Wtrial="<<std::endl;
	    std::cout<<"--------------:Wtrial n max_n_nodes=="<<Wtrial.node2elem().max_n_nodes()<<std::endl;

		// auto Waux=AuxFunctionSpacesBuild(pn);
		auto Waux=AuxFunctionSpacesBuild(pn,pn);
		std::cout<<"FIRST POST Waux="<<std::endl;

		auto W=FullSpaceBuild(Wtrial,Waux);
		std::cout<<" POST W="<<std::endl;
	    std::cout<<"--------------:W n max_n_nodes=="<<W.node2elem().max_n_nodes()<<std::endl;

		using W_type=decltype(W);
		auto W_ptr=std::make_shared<W_type>(W);

	    std::cout<<"--------------:W_ptr n max_n_nodes=="<<W_ptr->node2elem().max_n_nodes()<<std::endl;

		for(int i=0;i<n_levels;i++)
		{
		bisection.tracking_begin();
		bisection.uniform_refine(1);
		bisection.tracking_end();			
		}


		std::cout<<" BISECTION COMPLETED "<<std::endl;
 

	    std::cout<<"--------------:W_ptr n max_n_nodes=="<<W_ptr->node2elem().max_n_nodes()<<std::endl;



		auto sigma = MakeTrial<0>(W_ptr);
		auto u = MakeTrial<1>(W_ptr);

		auto tau = MakeTest<0>(W_ptr);
		auto v = MakeTest<1>(W_ptr);

		auto f1 = MakeFunction<0,ExactLinearElasticity<ManifoldDim>>(W_ptr);

		constexpr Real mu=1.0;
		constexpr Real lambda=1.0;
		constexpr Real alpha=1.0/Real(2.0*mu);
		constexpr Real beta=alpha*(-lambda/(Real(ManifoldDim)*lambda+ 2.0 * mu));

		constexpr auto halp=Constant<Scalar>(0.5);

		constexpr auto C_coeff1=Constant<Scalar>(1.0/Real(2.0*mu));

	    constexpr auto C_coeff2=Constant<Scalar>((1.0/Real(2.0*mu))*(-lambda/(Real(ManifoldDim)*lambda+ 2.0 * mu)));	
	    constexpr auto identity=Constant<Identity<ManifoldDim>>();					

	    auto Epsilon=NewOperator(halp*((GradientOperator()+Transpose(GradientOperator()))));//+Transpose(GradientOperator())));

		// auto C_inverse=NewOperator(C_coeff1*IdentityOperator() + C_coeff2 * IdentityOperator() * MatTraceOperator<Matrix<Real,2,2>>());
		auto C_inverse=NewOperator(C_coeff1*IdentityOperator() + C_coeff2*identity* MatTrace(IdentityOperator()) );

		// auto normal = MakeFunction<1,NormalFunction<ManifoldDim>,TraceOperator>(W_ptr);
		auto normal = MakeFunction<1,NormalFunction<ManifoldDim>>(W_ptr);
        // 

	    // 2D LSFEM POISSION
		auto bilinearform=
		// L2Inner(Div(sigma),Div(tau))+
		// L2Inner(C_inverse(sigma),C_inverse(tau))+
		// L2Inner(Epsilon(u),Epsilon(v))-
		// L2Inner(Epsilon(u),C_inverse(tau))-
		// L2Inner(C_inverse(sigma),Epsilon(v))

		// +surface_integral(0,Inner(Trace(sigma),normal),Inner(Trace(v),normal))
		+surface_integral(0,Trace(u),Trace(v))

		// L2Inner(MatTrace(sigma),MatTrace(tau))
		// L2Inner(sigma,halp*((Grad(v))))+
		// L2Inner(C_coeff2*identity* MatTrace(sigma),halp*(Transpose(Grad(v))))
		;

		auto linearform=
		L2Inner(-Div(tau),f1)
		// +
		// +surface_integral(0,Trace(u),(f1))
		;
        // decltype(MatTrace(IdentityOperator())) feef(5,5,5);
		std::cout<<" DEFINE BCS  "<<std::endl;




		auto bcs1=DirichletBC<1,FunctionZero<ManifoldDim>>(W_ptr,1);
		auto bcs2=DirichletBC<1,FunctionZero<ManifoldDim>>(W_ptr,2);
		auto bcs3=DirichletBC<1,FunctionZero<ManifoldDim>>(W_ptr,3);
		auto bcs4=DirichletBC<1,FunctionZero<ManifoldDim>>(W_ptr,4);
		auto bcs5=DirichletBC<1,FunctionZero<ManifoldDim>>(W_ptr,5);
		auto bcs6=DirichletBC<1,FunctionZero<ManifoldDim>>(W_ptr,6);


		auto context=create_context(bilinearform,linearform,bcs1,bcs2,bcs3,bcs4,bcs5,bcs6);

	    W_ptr->update();
	    
		SparseMatrix<Real> A;
		std::vector<Real> b;

		std::cout<<"ASSEMBLY"<<std::endl;
		std::cout<<"level---"<<level<<std::endl;
		context.assembly(A,b,level);

		// std::cout<<"APPLY BC "<<std::endl;
		// context.apply_bc(A,b);
		// std::vector<Real> x;
		// Integer max_iter=1;
  // 		std::ofstream os;
		// auto var_names=variables_names("stress","disp");

		// std::vector<Real> rhs;
  //       std::cout<<"START SOLVING PATCH MULTIGRID"<<std::endl;  
  //       gauss_seidel(x,A,b,1);         
		// std::cout<<"END SOLVING PATCH MULTIGRID"<<std::endl;
		
 	// 	std::string output_fileCOARSE ="LSFEM_ContactLinearElasticity"+ std::to_string(ManifoldDim) +
		// "D_RT" + std::to_string(Order1)+
		// "_P" + std::to_string(Order2)+"_outputCOARSE.vtk";

		// os.close();
  //       std::cout<<"GUARDA Qui"<<std::endl;
		// std::cout<<x.size()<<std::endl;
		// std::cout<<b.size()<<std::endl;
		// std::cout<<A.max_rows()<<std::endl;
		// std::cout<<A.max_cols()<<std::endl;
		// os.open(output_fileCOARSE.c_str());
		// write_wtk_isoparametric(os,W_ptr,x,var_names,level);


        


		// SparseMatrix<Real> AL;
		// std::vector<Real> bL;

		// Integer levelL=bisection.tracker().current_iterate()-1;

		// context.assembly(AL,bL,levelL);
	 //   // A.print_val();
		// std::cout<<"APPLY BC "<<std::endl;
		// context.apply_bc(AL,bL);


	 //   // A.print_val();
		// std::vector<Real> xL;


		// auto var_namesL=variables_names("stress","disp");


	 //     auto &dofsdofmap2=W_ptr->dofsdofmap();


    
	 //     auto& level_cumultive_n_dofs2=dofsdofmap2.level_cumultive_n_dofs();
	 //     std::cout<<"level_n_dofs_array="<<std::endl;
	 //     for(int i=0;i<level_cumultive_n_dofs2.size();i++)
	 //     {
	 //      // for(int j=0; j<level_n_dofs_array[i].size();j++)
	 //      std::cout<<level_cumultive_n_dofs2[i]<<" ";
	 //      std::cout<<std::endl;
	 //     }   

	 //    // FullFunctionSpaceInterpolation<decltype(W)> interp(W_ptr);
	 //    SparseMatrix<Real> InterpMat;
	 //    std::cout<<"level="<<level<<std::endl;
	 //    Integer level_finer=n_levels-1;
	 //    std::cout<<"level_finer="<<level_finer<<std::endl;
     
	 //    // interp.init(InterpMat,n_levels-1,n_levels);

	 //    FullFunctionSpaceLevelsInterpolation<W_type> levels_interp(W_ptr);
        
  //       std::vector<Integer> levels(n_levels+1-level);

  //       std::cout<<"n_levels=="<<n_levels <<std::endl;
  //       std::cout<<"level=="<<level <<std::endl;
  //       std::cout<<"n_levels-level=="<<n_levels-level <<std::endl;
  //       for(Integer i=0;i<n_levels+1-level;i++)
  //       {
  //       	levels[i]=i+level;
  //       	std::cout<<"levels[i]=="<<levels[i] <<std::endl;
  //       }

	 //    levels_interp.init(levels);


  //       std::cout<<"start multiply mat mat"<<std::endl;

       
  //       std::vector<SparseMatrix<Real>> A_levels(levels.size());

  //       // for(Integer el=0;el<mesh.n_elements();el++)
  //       // 	std::cout<<"el="<<el<<" lev="<<bisection.tracker().get_iterate(el)<<std::endl;
  //       // std::cout<<"i=="<<levels.size()-1<<std::endl;
  //       A_levels[levels.size()-1]=AL;

  //       std::cout<<A_levels[levels.size()-1].max_rows()<<std::endl;
  //       std::cout<<A_levels[levels.size()-1].max_cols()<<std::endl;

          
  //       std::cout<<" LEVEL INTERPS "<<std::endl;
  //       std::cout<<" LEVEL INTERPS "<<std::endl;


  //       for(Integer i=0;i<levels.size()-1;i++)
  //       {
  //       	auto& P=levels_interp.matrix(i);
  //       	std::cout<<P.max_rows()<<"  "<<P.max_cols()<<std::endl;
  //       }


  //       for(Integer i=levels.size()-2;i>=0;i--)
  //       { 
  //       	std::cout<<"i=="<<i<<std::endl;
  //       	auto& P=levels_interp.matrix(i);
  //           context.matrix_assembly(A_levels[i],levels[i]);
  //       	std::cout<<A_levels[i].max_rows()<<std::endl;
  //       	std::cout<<A_levels[i].max_cols()<<std::endl;
  //           context.build_boundary_info(levels[i]);
  //           context.apply_zero_bc_to_matrix(A_levels[i]);
  //       }


  //       std::vector<Real> solution;

  //       solution.resize(bL.size());

  //       for(Integer i=0;i<bL.size();i++)
  //       	solution[i]=0.0;


	 //    clock_t begin2 = clock();

	 //    Entity2Dofs<W_type,0> entity2dofs(W_ptr);
	 //    entity2dofs.build();
		// clock_t end2 = clock();
	 //    std::cout<<"TIME BUILD="<< double(end2 - begin2) / CLOCKS_PER_SEC <<std::endl;

  //        begin2 = clock();
	 //    auto& e2d=entity2dofs.get(levels);
	 //    end2 = clock();

	 //    std::cout<<"TIME E2D GET="<< double(end2 - begin2) / CLOCKS_PER_SEC <<std::endl;


	 //    clock_t begin1 = clock();
		// // auto& ordered_entities=entity2dofs.ordered_entities();	     
		// clock_t end1 = clock();
	 //    std::cout<<"TIME BUILD="<< double(end1 - begin1) / CLOCKS_PER_SEC <<std::endl;

  //         std::cout<<"END ORDERED ENTITY 2 DOFS"<<std::endl;
  //          // std::vector<Real> rhs;
           
  //        std::cout<<"START SOLVING PATCH MULTIGRID"<<std::endl;           
  //        patch_multigrid(solution,A_levels,bL,levels_interp,e2d,1,1,levels.size()-1,10,0.000000001);
  //        AL.multiply_and_add(rhs,-1.0,solution,bL);
  //        std::cout<<"residual="<<l2_norm(rhs)<<std::endl;
		// std::cout<<"END SOLVING PATCH MULTIGRID"<<std::endl;
		// std::string output_fileMULTIGRID ="LSFEM_ContactLinearElasticity"+ std::to_string(ManifoldDim) +
		// "D_RT" + std::to_string(Order1)+
		// "_P" + std::to_string(Order2)+"_outputMULTIGRID.vtk";

		// os.close();
		// os.open(output_fileMULTIGRID.c_str());
		// write_wtk_isoparametric(os,W_ptr,solution,var_namesL,levelL);
      
	}







	template<Integer ManifoldDim,Integer Order1,Integer Order2>
	void LSFEM_LinearElasticity(const Integer n, const Integer level, const Integer n_levels)
	{
	  // constexpr Integer ManifoldDim=2;
		constexpr Integer Dim=ManifoldDim;
		using MeshT=Mesh<Dim, ManifoldDim>;
		using Elem=typename MeshT::Elem;
		MeshT mesh;
		if(ManifoldDim==2)
		{
			std::cout<<"2D LSFEM_LinearElasticity";
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
			std::cout<<"3D LSFEM_LinearElasticity";
			generate_cube(mesh,n,n,n);
			mesh.build_dual_graph();
			mark_boundary(mesh);

	      // read_mesh("../data/tetrahedron.MFEM", mesh);
	      // assign_tags(mesh);
	      // Bisection<MeshT> bisection(mesh);
	      // bisection.uniform_refine(0);
	      // mesh.update_dual_graph();
		}
        // Integer n_levels=12;
	    clock_t begin=clock();
		clock_t end = clock();
		double elapsed_secs;

		// std::cout<<"mat==="<<mat3<<std::endl;

		

        std::cout<<"n2em init="<<std::endl;


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

	    Node2ElemMap<MeshT> n2em(mesh,bisection);
	    n2em.init();	      
	    // n2em.add_bisection(bisection);
	    // n2em.init();
	    std::cout<<"--------------:mesh n elements=="<<mesh.n_elements()<<std::endl;
	    std::cout<<"--------------:n2em n max_n_nodes=="<<n2em.max_n_nodes()<<std::endl;






		using AuxRT_n= FunctionSpace< MeshT, RT<Order1,ManifoldDim>>;
		using AuxP_n= FunctionSpace< MeshT, Lagrange<Order2,ManifoldDim>>;
		using AuxP_0= FunctionSpace< MeshT, Lagrange<0,ManifoldDim>>;
		using LSFEM= FunctionSpace< MeshT, RT<Order1,ManifoldDim>,Lagrange<Order2,ManifoldDim>>;

		AuxRT_n rtn(mesh,bisection,n2em);//csmc);
		AuxP_n pn(mesh,bisection,n2em);//csmc);
		AuxP_0 p0(mesh,bisection,n2em);//csmc);
		std::cout<<"FIRST PRE UPDATE="<<std::endl;
		LSFEM lsfem(mesh,bisection,n2em);//csmc);
		std::cout<<"FIRST POST UPDATE="<<std::endl;
	    std::cout<<"--------------:rtn n max_n_nodes=="<<rtn.node2elem().max_n_nodes()<<std::endl;
	    std::cout<<"--------------:pn n max_n_nodes=="<<pn.node2elem().max_n_nodes()<<std::endl;

	    auto Wtrial=MixedFunctionSpace(rtn,pn);
		std::cout<<"FIRST POST Wtrial="<<std::endl;
	    std::cout<<"--------------:Wtrial n max_n_nodes=="<<Wtrial.node2elem().max_n_nodes()<<std::endl;

		auto Waux=AuxFunctionSpacesBuild(pn);
		std::cout<<"FIRST POST Waux="<<std::endl;

		auto W=FullSpaceBuild(Wtrial,Waux);
		std::cout<<" POST W="<<std::endl;
	    std::cout<<"--------------:W n max_n_nodes=="<<W.node2elem().max_n_nodes()<<std::endl;

		using W_type=decltype(W);
		auto W_ptr=std::make_shared<W_type>(W);

	    std::cout<<"--------------:W_ptr n max_n_nodes=="<<W_ptr->node2elem().max_n_nodes()<<std::endl;

		for(int i=0;i<n_levels;i++)
		{
		bisection.tracking_begin();
		bisection.uniform_refine(1);
		bisection.tracking_end();			
		}


		std::cout<<" BISECTION COMPLETED "<<std::endl;
 

	    std::cout<<"--------------:W_ptr n max_n_nodes=="<<W_ptr->node2elem().max_n_nodes()<<std::endl;



		auto sigma = MakeTrial<0>(W_ptr);
		auto u = MakeTrial<1>(W_ptr);

		auto tau = MakeTest<0>(W_ptr);
		auto v = MakeTest<1>(W_ptr);

		auto f1 = MakeFunction<0,ExactLinearElasticity<ManifoldDim>>(W_ptr);

		constexpr Real mu=1.0;
		constexpr Real lambda=1.0;
		constexpr Real alpha=1.0/Real(2.0*mu);
		constexpr Real beta=alpha*(-lambda/(Real(ManifoldDim)*lambda+ 2.0 * mu));

		constexpr auto halp=Constant<Scalar>(0.5);

		constexpr auto C_coeff1=Constant<Scalar>(1.0/Real(2.0*mu));

	    constexpr auto C_coeff2=Constant<Scalar>((1.0/Real(2.0*mu))*(-lambda/(Real(ManifoldDim)*lambda+ 2.0 * mu)));	
	    constexpr auto identity=Constant<Identity<ManifoldDim>>();					

	    auto Epsilon=NewOperator(halp*((GradientOperator()+Transpose(GradientOperator()))));//+Transpose(GradientOperator())));

		// auto C_inverse=NewOperator(C_coeff1*IdentityOperator() + C_coeff2 * IdentityOperator() * MatTraceOperator<Matrix<Real,2,2>>());
		auto C_inverse=NewOperator(C_coeff1*IdentityOperator() + C_coeff2*identity* MatTrace(IdentityOperator()) );

		// auto normal = MakeFunction<1,NormalFunction<ManifoldDim>>(W_ptr);
        
        // auto normal =
        std::cout<<"C_coeff1="<<C_coeff1()<<std::endl;
        std::cout<<"C_coeff2="<<C_coeff2()<<std::endl;
        std::cout<<"identity="<<identity()<<std::endl;
		std::cout<<"RTN NDOFS="<<rtn.n_dofs()<<std::endl;
		std::cout<<"PN NDOFS="<<pn.n_dofs()<<std::endl;
		std::cout<<"NDOFS="<<W.spaces_ptr()->n_dofs()<<std::endl;

	    // 2D LSFEM POISSION
		auto bilinearform=
		L2Inner(Div(sigma),Div(tau))+
		L2Inner(C_inverse(sigma),C_inverse(tau))+
		L2Inner(Epsilon(u),Epsilon(v))-
		L2Inner(Epsilon(u),C_inverse(tau))-
		L2Inner(C_inverse(sigma),Epsilon(v))
		// L2Inner(MatTrace(sigma),MatTrace(tau))
		// L2Inner(sigma,halp*((Grad(v))))+
		// L2Inner(C_coeff2*identity* MatTrace(sigma),halp*(Transpose(Grad(v))))
		;

		auto linearform=
		L2Inner(-Div(tau),f1)
		;
        // decltype(MatTrace(IdentityOperator())) feef(5,5,5);
		std::cout<<" DEFINE BCS  "<<std::endl;




		auto bcs1=DirichletBC<1,FunctionZero<ManifoldDim>>(W_ptr,1);
		auto bcs2=DirichletBC<1,FunctionZero<ManifoldDim>>(W_ptr,2);
		auto bcs3=DirichletBC<1,FunctionZero<ManifoldDim>>(W_ptr,3);
		auto bcs4=DirichletBC<1,FunctionZero<ManifoldDim>>(W_ptr,4);
		auto bcs5=DirichletBC<1,FunctionZero<ManifoldDim>>(W_ptr,5);
		auto bcs6=DirichletBC<1,FunctionZero<ManifoldDim>>(W_ptr,6);

        




		// std::cout<<"CREATE CONTEXT"<<std::endl;
		// std::cout<<"coeff 2d "<<std::endl;
  //       std::cout<<"external "<<std::endl;
  //       // std::cout<<ReferenceShapeFunctionValue<Simplex<Dim,2>,IdentityOperator, RaviartThomasFE, 1>::reference_external_coeff<<std::endl;
  //       std::cout<<"internal "<<std::endl;

  //       // std::cout<<ReferenceShapeFunctionValue<Simplex<Dim,2>,IdentityOperator, RaviartThomasFE, 1>::reference_internal_coeff<<std::endl;

  //       std::cout<<"coeff 3d "<<std::endl;
  //       std::cout<<"external "<<std::endl;

  //       // std::cout<<ReferenceShapeFunctionValue<Simplex<Dim,3>,IdentityOperator, RaviartThomasFE, 1>::reference_external_coeff<<std::endl;
  //       std::cout<<"internal "<<std::endl;

        // std::cout<<ReferenceShapeFunctionValue<Simplex<Dim,3>,IdentityOperator, RaviartThomasFE, 1>::reference_internal_coeff<<std::endl;
        
        MapFromReference<IdentityOperator, Elem,RaviartThomasFE> map_RT;

        FiniteElem<Elem> FE(mesh);

        FE.init(2);
        typename ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, IdentityOperator, RaviartThomasFE, 1>::Output func;
  
        
        std::cout<<"apply3 1/3,1/3" <<std::endl;
        Vector<Real,2> point{1./3.,1./3.};
        // ReferenceShapeFunctionValue<Elem, IdentityOperator, RaviartThomasFE, 1>::apply3(func,point,FE,mesh);
        // std::cout<<func<<std::endl;
        std::cout<<"apply3 0,0" <<std::endl;
        point={0.,0.};
        // ReferenceShapeFunctionValue<Elem, IdentityOperator, RaviartThomasFE, 1>::apply3(func,point,FE,mesh);
        std::cout<<func<<std::endl;
        std::cout<<"apply3 1,0" <<std::endl;
        point={1.,0.};
        std::cout<<func<<std::endl;
        // ReferenceShapeFunctionValue<Elem, IdentityOperator, RaviartThomasFE, 1>::apply3(func,point,FE,mesh);
        point={0.,1.};
        std::cout<<"apply3 0,1" <<std::endl;
        // ReferenceShapeFunctionValue<Elem, IdentityOperator, RaviartThomasFE, 1>::apply3(func,point,FE,mesh);
        std::cout<<func<<std::endl;


		auto context=create_context(bilinearform,linearform,bcs1,bcs2,bcs3,bcs4,bcs5,bcs6);






	//     auto &dofsdofmap=W_ptr->dofsdofmap();
	// 	// auto &dofsdofmap=context.full_spaces_ptr()->dofsdofmap();
	// 	std::cout<<"<<<< UPDATE >>>>"<<std::endl;
	//     // dofsdofmap.update();
	//     std::cout<<"PRE UPDATE="<<std::endl;
	//     // context.full_spaces_ptr()->update();
	     W_ptr->update();
	    
	//     std::cout<<"POST UPDATE="<<std::endl;






		SparseMatrix<Real> A;
		std::vector<Real> b;

		std::cout<<"ASSEMBLY"<<std::endl;
		// Integer level=2;//bisection.tracker().current_iterate();
		// Integer level=4;
		std::cout<<"level---"<<level<<std::endl;
		context.assembly(A,b,level);


		// A.print_val();
	//    // A.print_val();
		std::cout<<"APPLY BC "<<std::endl;
		context.apply_bc(A,b);


	//    // A.print_val();
		std::vector<Real> x;
		Integer max_iter=1;
  		std::ofstream os;
		auto var_names=variables_names("stress","disp");









std::vector<Real> rhs;
         std::cout<<"START SOLVING PATCH MULTIGRID"<<std::endl;  
         gauss_seidel(x,A,b,1);         
		std::cout<<"END SOLVING PATCH MULTIGRID"<<std::endl;
		
  //        A.multiply_and_add(rhs,-1.0,x,b);
  //        std::cout<<"residual="<<l2_norm(rhs)<<std::endl;


		std::string output_fileCOARSE ="LSFEM_Elasticity"+ std::to_string(ManifoldDim) +
		"D_RT" + std::to_string(Order1)+
		"_P" + std::to_string(Order2)+"_outputCOARSE.vtk";

		os.close();
        std::cout<<"GUARDA Qui"<<std::endl;
		std::cout<<x.size()<<std::endl;
		std::cout<<b.size()<<std::endl;
		std::cout<<A.max_rows()<<std::endl;
		std::cout<<A.max_cols()<<std::endl;
		os.open(output_fileCOARSE.c_str());
		write_wtk_isoparametric(os,W_ptr,x,var_names,level);

		// std::cout<<"START SOLVING"<<std::endl;
		// // gauss_seidel(x,A,b,max_iter);
		// std::cout<<"END SOLVING"<<std::endl;

	    // std::cout<<"---------------------------x"<<std::endl;





        


		SparseMatrix<Real> AL;
		std::vector<Real> bL;

		Integer levelL=bisection.tracker().current_iterate()-1;

		context.assembly(AL,bL,levelL);
	   // A.print_val();
		std::cout<<"APPLY BC "<<std::endl;
		context.apply_bc(AL,bL);


	   // A.print_val();
		std::vector<Real> xL;


		auto var_namesL=variables_names("stress","disp");


	     auto &dofsdofmap2=W_ptr->dofsdofmap();


    
	     auto& level_cumultive_n_dofs2=dofsdofmap2.level_cumultive_n_dofs();
	     std::cout<<"level_n_dofs_array="<<std::endl;
	     for(int i=0;i<level_cumultive_n_dofs2.size();i++)
	     {
	      // for(int j=0; j<level_n_dofs_array[i].size();j++)
	      std::cout<<level_cumultive_n_dofs2[i]<<" ";
	      std::cout<<std::endl;
	     }   

	    // FullFunctionSpaceInterpolation<decltype(W)> interp(W_ptr);
	    SparseMatrix<Real> InterpMat;
	    std::cout<<"level="<<level<<std::endl;
	    Integer level_finer=n_levels-1;
	    std::cout<<"level_finer="<<level_finer<<std::endl;
     
	    // interp.init(InterpMat,n_levels-1,n_levels);

	    FullFunctionSpaceLevelsInterpolation<W_type> levels_interp(W_ptr);
        
        std::vector<Integer> levels(n_levels+1-level);

        std::cout<<"n_levels=="<<n_levels <<std::endl;
        std::cout<<"level=="<<level <<std::endl;
        std::cout<<"n_levels-level=="<<n_levels-level <<std::endl;
        for(Integer i=0;i<n_levels+1-level;i++)
        {
        	levels[i]=i+level;
        	std::cout<<"levels[i]=="<<levels[i] <<std::endl;
        }

	    levels_interp.init(levels);


        std::cout<<"start multiply mat mat"<<std::endl;

       
        std::vector<SparseMatrix<Real>> A_levels(levels.size());

        // for(Integer el=0;el<mesh.n_elements();el++)
        // 	std::cout<<"el="<<el<<" lev="<<bisection.tracker().get_iterate(el)<<std::endl;
        // std::cout<<"i=="<<levels.size()-1<<std::endl;
        A_levels[levels.size()-1]=AL;

        std::cout<<A_levels[levels.size()-1].max_rows()<<std::endl;
        std::cout<<A_levels[levels.size()-1].max_cols()<<std::endl;

          
        std::cout<<" LEVEL INTERPS "<<std::endl;
        std::cout<<" LEVEL INTERPS "<<std::endl;


        for(Integer i=0;i<levels.size()-1;i++)
        {
        	auto& P=levels_interp.matrix(i);
        	std::cout<<P.max_rows()<<"  "<<P.max_cols()<<std::endl;
        }


        for(Integer i=levels.size()-2;i>=0;i--)
        { 
        	std::cout<<"i=="<<i<<std::endl;
        	auto& P=levels_interp.matrix(i);
        	// auto AP=A_levels[i+1].multiply(levels_interp.matrix(i));
        	// std::cout<<P.max_rows()<<std::endl;
        	// std::cout<<P.max_cols()<<std::endl;
        	// std::cout<<AP.max_rows()<<std::endl;
        	// std::cout<<AP.max_cols()<<std::endl;
         //    A_levels[i]=levels_interp.matrix(i).transpose_and_multiply(AP);

            context.matrix_assembly(A_levels[i],levels[i]);
        	std::cout<<A_levels[i].max_rows()<<std::endl;
        	std::cout<<A_levels[i].max_cols()<<std::endl;
            context.build_boundary_info(levels[i]);
            context.apply_zero_bc_to_matrix(A_levels[i]);
        }


        std::vector<Real> solution;

        solution.resize(bL.size());

        for(Integer i=0;i<bL.size();i++)
        	solution[i]=0.0;


	    clock_t begin2 = clock();

	    Entity2Dofs<W_type,0> entity2dofs(W_ptr);
	    entity2dofs.build();
		clock_t end2 = clock();
	    std::cout<<"TIME BUILD="<< double(end2 - begin2) / CLOCKS_PER_SEC <<std::endl;

         begin2 = clock();
	    auto& e2d=entity2dofs.get(levels);
	    end2 = clock();

	    std::cout<<"TIME E2D GET="<< double(end2 - begin2) / CLOCKS_PER_SEC <<std::endl;


	    clock_t begin1 = clock();
		// auto& ordered_entities=entity2dofs.ordered_entities();	     
		clock_t end1 = clock();
	    std::cout<<"TIME BUILD="<< double(end1 - begin1) / CLOCKS_PER_SEC <<std::endl;

          std::cout<<"END ORDERED ENTITY 2 DOFS"<<std::endl;
           // std::vector<Real> rhs;
           
         std::cout<<"START SOLVING PATCH MULTIGRID"<<std::endl;           
         patch_multigrid(solution,A_levels,bL,levels_interp,e2d,1,1,levels.size()-1,10,0.000000001);
         AL.multiply_and_add(rhs,-1.0,solution,bL);
         std::cout<<"residual="<<l2_norm(rhs)<<std::endl;
		std::cout<<"END SOLVING PATCH MULTIGRID"<<std::endl;
		std::string output_fileMULTIGRID ="LSFEM_Elasticity"+ std::to_string(ManifoldDim) +
		"D_RT" + std::to_string(Order1)+
		"_P" + std::to_string(Order2)+"_outputMULTIGRID.vtk";

		os.close();
		os.open(output_fileMULTIGRID.c_str());
		write_wtk_isoparametric(os,W_ptr,solution,var_namesL,levelL);

		// Matrix<Real,4,2> points{0,0,1.,0,0,1.,1./3,1./3};
		// Array<Real,8> alpha2{1,1,1,1,1,1,1,1};

		// FQPValues<Matrix<Real,2,1>,4,8> func_values;
		// FQPValues<Matrix<Real,1,1>,4,8> div_values;
		// MapFromReference<IdentityOperator,Elem,RaviartThomasFE> map_rt;
		// MapFromReference<DivergenceOperator,Elem,RaviartThomasFE> map_div_rt;

		// FE.init(2);
		// map_rt.init(FE);
		// map_div_rt.init(FE);

		// // std::cout<<"IDENTITY dynamic"<<std::endl;   
  // //       DynamicShapeFunctionValue<ElementFunctionSpace<Elem,RaviartThomasFE,1,1,1>,IdentityOperator>::apply4(func_values,points,FE,map_rt(),alpha2); 

  // //       ShapeFunction<Simplex<2,2>,BaseFunctionSpace<RaviartThomasFE,1,1,1>,IdentityOperator,GaussPoints< Simplex<Dim,2> , 3>> shape;
  // //       std::cout<<"IDENTITY reference values"<<std::endl;

  // //       shape.init_map(map_rt);
  // //       shape.init(alpha2,points,FE);
  // //       std::cout<<"after init"<<std::endl;
  // //       std::cout<<shape.eval()<<std::endl;


		// std::cout<<"DIVERGENCE dynamic"<<std::endl;   
  //       DynamicShapeFunctionValue<ElementFunctionSpace<Elem,RaviartThomasFE,1,1,1>,DivergenceOperator>::apply4(div_values,points,FE,map_div_rt(),alpha2); 

  //       ShapeFunction<Simplex<2,2>,BaseFunctionSpace<RaviartThomasFE,1,1,1>,DivergenceOperator,GaussPoints< Simplex<Dim,2> , 3>> shape_div;
  //       std::cout<<"DIVERGENCE reference values"<<std::endl;

  //       shape_div.init_map(map_div_rt);
  //       shape_div.init(alpha2,points,FE);
  //       std::cout<<"after init"<<std::endl;
  //       std::cout<<shape_div.eval()<<std::endl;
      


      DofAux<ElementFunctionSpace<Simplex<Dim,ManifoldDim>,RaviartThomasFE,1,1,1>> dof ;

      FiniteElem<Elem> C_FE(mesh);
      FiniteElem<Elem> F_FE(mesh);

      Matrix<Real,8,8> mat88;

      C_FE.init(0);
      F_FE.init(2);




      dof.compute(mat88,C_FE,F_FE);


	    std::vector<std::vector<Real>> x_interp_vec(levels.size());

        std::cout<<"multiply level == "<<0<<std::endl;
	    x_interp_vec[0]=x;

        for(Integer i=1;i<levels.size();i++)
        {
            // std::cout<<"pre ultiply level == "<<i<<std::endl;
        	x_interp_vec[i]=levels_interp.matrix(i-1).multiply(x_interp_vec[i-1]);
        	// std::cout<<"after multiply level == "<<i<<std::endl;
			std::string output_file_interp ="LSFEM_Elasticity"+ std::to_string(ManifoldDim) +
			"D_RT" + std::to_string(Order1)+
			"_P" + std::to_string(Order2)+"_outputINTERPLEVEL"+ std::to_string(levels[i])+"-"+
			std::to_string(levels[i+1])+
			+".vtk";
			os.close();
			os.open(output_file_interp.c_str());
			// write_wtk_isoparametric(os,*context.full_spaces_ptr(),xL,var_namesL,level);
			write_wtk_isoparametric(os,W_ptr,x_interp_vec[i],var_namesL,levels[i]);
		    os.close();
        }



	}



























































	template<Integer ManifoldDim,Integer Order1,Integer Order2>
	void LSFEM_Poisson(const Integer n, const Integer level, const Integer n_levels)
	{
	  // constexpr Integer ManifoldDim=2;
		constexpr Integer Dim=ManifoldDim;
		using MeshT=Mesh<Dim, ManifoldDim>;
		using Elem=typename MeshT::Elem;
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
        // Integer n_levels=12;
	    clock_t begin=clock();
		clock_t end = clock();
		double elapsed_secs;

		// std::cout<<"mat==="<<mat3<<std::endl;

		

        std::cout<<"n2em init="<<std::endl;


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
	    
	  
	    Node2ElemMap<MeshT> n2em(mesh,bisection);
	    n2em.init();	      
	    // n2em.add_bisection(bisection);
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
		using AuxP_0= FunctionSpace< MeshT, Lagrange<0,1>>;
		using LSFEM= FunctionSpace< MeshT, RT<Order1,1>,Lagrange<Order2,1>>;

	    // ConnectivitySimpliacialMapCollection<MeshT> csmc(mesh,n2em,bisection);
	    // csmc.template init<RT<Order1,1>,Lagrange<Order2,1>>();

	    // csmc.describe2(std::cout);
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
AuxP_0 p0(mesh,bisection,n2em);//csmc);
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
std::cout<<"FIRST POST Wtrial="<<std::endl;

		auto Waux=AuxFunctionSpacesBuild(pn);
std::cout<<"FIRST POST Waux="<<std::endl;

		auto W=FullSpaceBuild(Wtrial,Waux);
std::cout<<" POST W="<<std::endl;

		using W_type=decltype(W);
		auto W_ptr=std::make_shared<W_type>(W);


		for(int i=0;i<n_levels;i++)
		{
		bisection.tracking_begin();
		bisection.uniform_refine(1);
		bisection.tracking_end();			
		}


		std::cout<<" BISECTION COMPLETED "<<std::endl;
 




		auto u0 = MakeTrial<0>(W_ptr);
		auto u1 = MakeTrial<1>(W_ptr);

		auto v0 = MakeTest<0>(W_ptr);
		auto v1 = MakeTest<1>(W_ptr);

		auto f1 = MakeFunction<0,ExactPoisson<ManifoldDim>>(W_ptr);
		// auto normal = MakeFunction<1,NormalFunction<ManifoldDim>>(W_ptr);
        
        // auto normal =
		std::cout<<"NDOFS="<<W.spaces_ptr()->n_dofs()<<std::endl;

	  // 2D LSFEM POISSION
		auto bilinearform=L2Inner(Div(u0),Div(v0))+
		L2Inner((u0),(v0))+
		L2Inner(Grad(u1),Grad(v1))-
		L2Inner(Grad(u1),(v0))-
		L2Inner((u0),Grad(v1))
		;
	  // auto bilinearform=//L2Inner(Div(u0),Div(v0))+
	  //                   L2Inner((u0),(v0));   

		auto linearform=
		L2Inner(-Div(v0),f1)
		// +surface_integral(0,normal,v0)
		;

		std::cout<<" DEFINE BCS  "<<std::endl;




		auto bcs1=DirichletBC<1,FunctionZero<ManifoldDim>>(W_ptr,1);
		auto bcs2=DirichletBC<1,FunctionZero<ManifoldDim>>(W_ptr,2);
		auto bcs3=DirichletBC<1,FunctionZero<ManifoldDim>>(W_ptr,3);
		auto bcs4=DirichletBC<1,FunctionZero<ManifoldDim>>(W_ptr,4);
		auto bcs5=DirichletBC<1,FunctionZero<ManifoldDim>>(W_ptr,5);
		auto bcs6=DirichletBC<1,FunctionZero<ManifoldDim>>(W_ptr,6);

        
		// using P1_2= FunctionSpace< MeshT, Lagrange<1,2>>;
		// P1_2 p12(mesh,bisection,n2em);
		// auto Xtrial=MixedFunctionSpace(p12);
		// auto Xaux=AuxFunctionSpacesBuild(p0);
		// auto X=FullSpaceBuild(Xtrial,Xaux);
		// auto X_ptr=std::make_shared<decltype(X)>(X);
		// auto x0 = MakeTest<0>(X_ptr);
		// auto normal = MakeFunction<0,NormalFunction<ManifoldDim>>(X_ptr);


  //       auto tryform=surface_integral(0,normal,Trace(x0));
  //       auto generaltryform=general_form(tryform);






		std::cout<<"CREATE CONTEXT"<<std::endl;
        auto genform=general_form(bilinearform);

		std::cout<<"bilinearform->"<<genform.spaces_ptr()->node2elem().max_n_nodes()<<std::endl;
		std::cout<<"bilinearform->"<<genform.spaces_ptr()->spaces_ptr()->node2elem().max_n_nodes()<<std::endl;


		auto context=create_context(bilinearform,linearform,bcs1,bcs2,bcs3,bcs4,bcs5,bcs6);





	    auto &dofsdofmap=W_ptr->dofsdofmap();
		// auto &dofsdofmap=context.full_spaces_ptr()->dofsdofmap();
		std::cout<<"<<<< UPDATE >>>>"<<std::endl;
	    // dofsdofmap.update();
	    std::cout<<"PRE UPDATE="<<std::endl;
	    // context.full_spaces_ptr()->update();
	     W_ptr->update();
	    
	    // auto& level_cumultive_n_dofs=dofsdofmap.level_cumultive_n_dofs();
	    //  std::cout<<"level_n_dofs_array="<<std::endl;
	    //  for(int i=0;i<level_cumultive_n_dofs.size();i++)
	    //  {
	    //   // for(int j=0; j<level_n_dofs_array[i].size();j++)
	    //     std::cout<<level_cumultive_n_dofs[i]<<" ";
	    //   std::cout<<std::endl;
	    //  }   
	    std::cout<<"POST UPDATE="<<std::endl;






		SparseMatrix<Real> A;
		std::vector<Real> b;

		std::cout<<"ASSEMBLY"<<std::endl;
		// Integer level=2;//bisection.tracker().current_iterate();
		// Integer level=4;
		std::cout<<"level---"<<level<<std::endl;
		context.assembly(A,b,level);
	   // A.print_val();
		std::cout<<"APPLY BC "<<std::endl;
		context.apply_bc(A,b);


	   // A.print_val();
		std::vector<Real> x;
		Integer max_iter=1;


		std::cout<<"START SOLVING"<<std::endl;
		// gauss_seidel(x,A,b,max_iter);
		std::cout<<"END SOLVING"<<std::endl;

	    std::cout<<"---------------------------x"<<std::endl;
	    ////// MODIFY X BEFORE INTERPOLATIONA
	    // for(int i=0;i<x.size();i++)
	    	// x[i]=1;


	    // for(int i=0;i<x.size();i++)
	    // 	std::cout<<x[i]<<std::endl;





		std::ofstream os;
		auto var_names=variables_names("stress","disp");
		// std::string output_file ="LSFEM_Poisson"+ std::to_string(ManifoldDim) +
		// "D_RT" + std::to_string(Order1)+
		// "_P" + std::to_string(Order2)+"_output.vtk";

		// os.close();
		// os.open(output_file.c_str());
		// // write_wtk_isoparametric(os,*context.full_spaces_ptr(),x,var_names,level);
		// std::cout<<" write to file START"<<std::endl;
		// write_wtk_isoparametric(os,W_ptr,x,var_names,level);
		// std::cout<<" write to file END"<<std::endl;

	 //    os.close();




     //    std::cout<<std::endl;
	    // std::cout<<"----------b"<<std::endl;
	    // for(int i=0;i<b.size();i++)
	    // 	std::cout<<b[i]<<std::endl;



        


		SparseMatrix<Real> AL;
		std::vector<Real> bL;

		Integer levelL=bisection.tracker().current_iterate()-1;

		context.assembly(AL,bL,levelL);
	   // A.print_val();
		std::cout<<"APPLY BC "<<std::endl;
		context.apply_bc(AL,bL);


	   // A.print_val();
		std::vector<Real> xL;


		// std::cout<<"START SOLVING"<<std::endl;
		// gauss_seidel(xL,AL,bL,max_iter);
		// std::cout<<"END SOLVING"<<std::endl;
		auto var_namesL=variables_names("stress","disp");
		// std::string output_fileL ="LSFEM_Poisson"+ std::to_string(ManifoldDim) +
		// "D_RT" + std::to_string(Order1)+
		// "_P" + std::to_string(Order2)+"_outputFINE.vtk";

		// os.close();
		// os.open(output_fileL.c_str());
		// // write_wtk_isoparametric(os,*context.full_spaces_ptr(),xL,var_namesL,level);
		// write_wtk_isoparametric(os,W_ptr,xL,var_namesL,levelL);

	 //    os.close();


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

	    // FullFunctionSpaceInterpolation<decltype(W)> interp(W_ptr);
	    SparseMatrix<Real> InterpMat;
	    std::cout<<"level="<<level<<std::endl;
	    Integer level_finer=n_levels-1;
	    std::cout<<"level_finer="<<level_finer<<std::endl;
     
	    // interp.init(InterpMat,n_levels-1,n_levels);

	    FullFunctionSpaceLevelsInterpolation<W_type> levels_interp(W_ptr);
        
        std::vector<Integer> levels(n_levels+1-level);

        std::cout<<"n_levels=="<<n_levels <<std::endl;
        std::cout<<"level=="<<level <<std::endl;
        std::cout<<"n_levels-level=="<<n_levels-level <<std::endl;
        for(Integer i=0;i<n_levels+1-level;i++)
        {
        	levels[i]=i+level;
        	std::cout<<"levels[i]=="<<levels[i] <<std::endl;
        }

	    levels_interp.init(levels);

	  //   std::vector<std::vector<Real>> x_interp_vec(levels.size());

   //      std::cout<<"multiply level == "<<0<<std::endl;
	  //   x_interp_vec[0]=x;

   //      for(Integer i=1;i<levels.size();i++)
   //      {
   //          // std::cout<<"pre ultiply level == "<<i<<std::endl;
   //      	x_interp_vec[i]=levels_interp.matrix(i-1).multiply(x_interp_vec[i-1]);
   //      	// std::cout<<"after multiply level == "<<i<<std::endl;
			// std::string output_file_interp ="LSFEM_Poisson"+ std::to_string(ManifoldDim) +
			// "D_RT" + std::to_string(Order1)+
			// "_P" + std::to_string(Order2)+"_outputINTERPLEVEL"+ std::to_string(levels[i])+"-"+
			// std::to_string(levels[i+1])+
			// +".vtk";
			// os.close();
			// os.open(output_file_interp.c_str());
			// // write_wtk_isoparametric(os,*context.full_spaces_ptr(),xL,var_namesL,level);
			// write_wtk_isoparametric(os,W_ptr,x_interp_vec[i],var_namesL,levels[i]);
		 //    os.close();
   //      }

	    // InterpMat.print_val();


	    // std::cout<<"before x_interp"<<std::endl;
	    // auto x_interp=InterpMat.multiply(x);
	    
//     std::vector<Real> x_interp
// {-6.38588e-16,-0.8112,1.52498e-16,-7.45169e-16,-0.8112,0.4056,3.26509e-16,0.4056,-0.8112,-0.8112,-2.67354e-16,-0.8112,-0.8112,-0.8112,-9.00613e-19,-0.8112,0.4056,1.52498e-16,0.4056,-2.72351e-16,-2.72351e-16,0.4056,0.4056,3.26509e-16,0.4056,-4.35486e-16,-4.35486e-16,0.4056,0,0,0,0,0.5408,0,0,0,0,-0.2704,-0.2704,-0.2704,-0.2704};
//         x_interp[181]=0.1933;

    	// std::cout<<"x_interp"<<std::endl;
	    // for(int i=0;i<x_interp.size();i++)
	    // 	std::cout<<x_interp[i]<<std::endl;


	    // x_interp[180]=-0.1933;
		// std::string output_file_interp ="LSFEM_Poisson"+ std::to_string(ManifoldDim) +
		// "D_RT" + std::to_string(Order1)+
		// "_P" + std::to_string(Order2)+"_outputINTERP.vtk";

		// os.close();
		// os.open(output_file_interp.c_str());
		// // write_wtk_isoparametric(os,*context.full_spaces_ptr(),xL,var_namesL,level);
		// write_wtk_isoparametric(os,W_ptr,x_interp,var_namesL,level_finer);

	 //    os.close();


        std::cout<<"start multiply mat mat"<<std::endl;

       
        std::vector<SparseMatrix<Real>> A_levels(levels.size());


        // for(Integer el=0;el<mesh.n_elements();el++)
        	// std::cout<<"el="<<el<<" lev="<<bisection.tracker().get_iterate(el)<<std::endl;
        // std::cout<<"i=="<<levels.size()-1<<std::endl;
        A_levels[levels.size()-1]=AL;

        std::cout<<A_levels[levels.size()-1].max_rows()<<std::endl;
        std::cout<<A_levels[levels.size()-1].max_cols()<<std::endl;

          
        std::cout<<" LEVEL INTERPS "<<std::endl;
        std::cout<<" LEVEL INTERPS "<<std::endl;


        for(Integer i=0;i<levels.size()-1;i++)
        {
        	auto& P=levels_interp.matrix(i);
        	std::cout<<P.max_rows()<<"  "<<P.max_cols()<<std::endl;
        }


        for(Integer i=levels.size()-2;i>=0;i--)
        { 
        	std::cout<<"i=="<<i<<std::endl;
        	auto& P=levels_interp.matrix(i);
        	// auto AP=A_levels[i+1].multiply(levels_interp.matrix(i));
        	// std::cout<<P.max_rows()<<std::endl;
        	// std::cout<<P.max_cols()<<std::endl;
        	// std::cout<<AP.max_rows()<<std::endl;
        	// std::cout<<AP.max_cols()<<std::endl;
         //    A_levels[i]=levels_interp.matrix(i).transpose_and_multiply(AP);

            context.matrix_assembly(A_levels[i],levels[i]);
        	std::cout<<A_levels[i].max_rows()<<std::endl;
        	std::cout<<A_levels[i].max_cols()<<std::endl;
            context.build_boundary_info(levels[i]);
            context.apply_zero_bc_to_matrix(A_levels[i]);
            // A_levels[i].print_val();
        }

        // for(Integer i=0;i<A_levels.size();i++)
        // {
        // 	A_levels[i].print_val();
        // }

        // for(Integer i=0;i<A_levels.size()-1;i++)
        // {
        // 	levels_interp.matrix(i).print_val();
        // }

        std::vector<Real> solution;

        solution.resize(bL.size());

        for(Integer i=0;i<bL.size();i++)
        	solution[i]=0.0;

  //       multigrid(solution,A_levels,bL,levels_interp,3,3,levels.size()-1,10,0.0000000001);

		// std::string output_fileMULTIGRID ="LSFEM_Poisson"+ std::to_string(ManifoldDim) +
		// "D_RT" + std::to_string(Order1)+
		// "_P" + std::to_string(Order2)+"_outputMULTIGRID.vtk";

		// os.close();
		// os.open(output_fileMULTIGRID.c_str());
		// // write_wtk_isoparametric(os,*context.full_spaces_ptr(),x,var_names,level);
		// write_wtk_isoparametric(os,W_ptr,solution,var_names,levels[levels.size()-1]);

	 //    os.close();









		// using AuxRT1= FunctionSpace< MeshT, RT<0,2>>;
		// using AuxP2= FunctionSpace< MeshT, Lagrange<1,1>>;

		// AuxRT1 rt1(mesh,bisection,n2em);//csmc);
		// AuxP2 p2(mesh,bisection,n2em);//csmc);

	 //    auto Wtrial2try=(MixedFunctionSpace(rt1,p2));
	 //    auto Wtrial2=FullSpaceBuild(Wtrial2try);
	 //    using W_type2=decltype(Wtrial2);
	 //    auto W_ptr2=std::make_shared<W_type2>(Wtrial2);


	    clock_t begin2 = clock();

	    Entity2Dofs<W_type,0> entity2dofs(W_ptr);
	    entity2dofs.build();
		clock_t end2 = clock();
	    std::cout<<"TIME BUILD="<< double(end2 - begin2) / CLOCKS_PER_SEC <<std::endl;

         begin2 = clock();
	    auto& e2d=entity2dofs.get(levels);
	    end2 = clock();

	    std::cout<<"TIME E2D GET="<< double(end2 - begin2) / CLOCKS_PER_SEC <<std::endl;
	    // entity2dofs.print();

	   //  auto& e0=entity2dofs.get(0);

	    
	   //  auto& e1=entity2dofs.get(1);
	   //    for(Integer i=0;i<e0.size();i++)
	   //    {
	   //    	for(Integer j=0;j<e0[i].size();j++)
	   //        std::cout<<e0[i][j]<<" ";
	   //       std::cout<< std::endl;
	   //    }              
	   //      std::cout<< std::endl;
	   //      std::cout<< std::endl;
	   //    for(Integer i=0;i<e1.size();i++)
	   //    {
	   //    	for(Integer j=0;j<e1[i].size();j++)
	   //        std::cout<<e1[i][j]<<" ";
	   //       std::cout<< std::endl;
	   //    }  

    //       std::cout<<"START ENTITY 2 DOFS"<<std::endl;
	   //    for(Integer k=0;k<e2d.size();k++)
	   //    {
	   //    	 std::cout<<"level = "<<k <<std::endl;
		  //     for(Integer i=0;i<e2d[k].size();i++)
		  //     {
		  //     	for(Integer j=0;j<e2d[k][i].size();j++)
		  //         std::cout<<e2d[k][i][j]<<" ";
		  //        std::cout<< std::endl;
		  //     } 
		  //      // std::cout<< std::endl;
	 	 // }
	 	 // std::cout<<"END ENTITY 2 DOFS"<<std::endl;
	




	     // densemat.init(mat4);
	     // densemat.cholesky_matrix(mat4_tmp);


	     // std::cout<<"cholesky mat"<<std::endl;

	     // densemat.cholesky(sol4,vec4);
      //    for(Integer i=0;i<4;i++)
      //    std::cout<<sol4[i]<<std::endl;
         // AL.print_val();


         // patch_gauss_seidel(solutionL,AL,bL,e1,max_iter);


	     

 	    // InterpMat.print_val();
	    // AL.print_val();

	 //    auto AP=AL.multiply(InterpMat);

	 //    // AP.print_val();

	 //    auto PAP=InterpMat.transpose_and_multiply(AP);

	 //    // PAP.print_val();


		// SparseMatrix<Real> AProva;
		// std::vector<Real> bprova;

		// context.assembly(AProva,bprova,levelL-1);
	 //   // A.print_val();
		// std::cout<<"APPLY BC "<<std::endl;
		// context.apply_bc(PAP,bprova);


		// std::vector<Real> xprova;


		// std::cout<<"START SOLVING"<<std::endl;
		// gauss_seidel(xprova,PAP,bprova,max_iter);
		// std::cout<<"END SOLVING"<<std::endl;
		// std::string output_fileprova ="LSFEM_Poisson"+ std::to_string(ManifoldDim) +
		// "D_RT" + std::to_string(Order1)+
		// "_P" + std::to_string(Order2)+"_outputPROVA.vtk";

		// os.close();
		// os.open(output_fileprova.c_str());
		// // write_wtk_isoparametric(os,*context.full_spaces_ptr(),xL,var_namesL,level);
		// write_wtk_isoparametric(os,W_ptr,xprova,var_names,levelL-1);

	 //    os.close();




	    // SparseMatrix<Real> mat;

	    // mat.init(5,5,5);

	    // mat.equal(10,0,0);
	    // mat.equal(10,0,3);
	    // mat.equal(10,0,4);

	    // mat.equal(10,1,0);
	    // mat.equal(10,1,4);

	    // mat.equal(10,2,0);
	    // mat.equal(10,2,1);
	    // mat.equal(10,2,3);
	    // mat.equal(10,2,4);

	    // mat.equal(10,3,2);
	    // mat.equal(10,3,3);
	    // mat.equal(10,3,4);

	    // mat.equal(10,4,0);
	    // mat.equal(10,4,1);
	    // mat.equal(10,4,2);
	    // mat.equal(10,4,3);
	    // mat.equal(10,4,4);

	    // auto& m1=mat.cols_idx();
	    // auto& m2=mat.rows_idx();

	    // std::vector<Real> vec6{1.2,2.3,3.4,4.5,5.6};

	    // auto vec7=mat.transpose_and_multiply(vec6);

	    // std::cout<<"vec6"<<std::endl;
	    // for(Integer i=0;i<vec7.size();i++)
	    // 	std::cout<<vec6[i]<<std::endl;

	    // std::cout<<"vec7"<<std::endl;
	    // for(Integer i=0;i<vec7.size();i++)
	    // 	std::cout<<vec7[i]<<std::endl;
	    // mat.print_val();












	    // AP.print_val();






   //  std::cout<<"LagrangeDofs1="<<std::endl;
   //  std::cout<<LagrangeDofs<ElemFunctionSpace<Simplex<3,3>,Lagrange1<1>> >::value()<<std::endl;
   // std::cout<<"LagrangeDofs1="<<std::endl;

   //  std::cout<<LagrangeDofs<ElemFunctionSpace<Simplex<4,3>,Lagrange1<1>> >::value()<<std::endl;
   //  std::cout<<"LagrangeDofs1="<<std::endl;
   //   std::cout<<LagrangeDofs<ElemFunctionSpace<Simplex<4,4>,Lagrange1<1>> >::value()<<std::endl;

 


   //  std::cout<<"LagrangeDofs2="<<std::endl;
   //  std::cout<<LagrangeDofs<ElemFunctionSpace<Simplex<3,3>,Lagrange2<1>> >::value()<<std::endl;
   //  std::cout<<"LagrangeDofs2="<<std::endl;
   //  std::cout<<LagrangeDofs<ElemFunctionSpace<Simplex<4,3>,Lagrange2<1>> >::value()<<std::endl;
   //  std::cout<<"LagrangeDofs2="<<std::endl;
   //  std::cout<<LagrangeDofs<ElemFunctionSpace<Simplex<4,4>,Lagrange2<1>> >::value()<<std::endl;
   //  Dof<ElemFunctionSpace<Simplex<Dim,ManifoldDim>,Lagrange1<1>>> dof_p1(mesh);
   //  dof_p1.compute(mesh.elem(0),mesh.elem(mesh.elem(0).children[0]));
   //  // dofp1.compute(mesh.elem(0),mesh.elem(3));
   //  // dofp1.compute(mesh.elem(1),mesh.elem(4));
   //  // dofp1.compute(mesh.elem(1),mesh.elem(5));

   //  DofAux<ElemFunctionSpace<Simplex<Dim,ManifoldDim>,Lagrange1<1>>> dofp1;
   //  DofAux<ElemFunctionSpace<Simplex<Dim,ManifoldDim>,Lagrange2<1>>> dofp2;
   //  DofAux<ElemFunctionSpace<Simplex<Dim,ManifoldDim>,RT0<1>>> dofrt0;
   //  DofAux<ElemFunctionSpace<Simplex<Dim,ManifoldDim>,RT1<1>>> dofrt1;
   //  FiniteElem<Elem> C_FE(mesh);
   //  FiniteElem<Elem> F_FE(mesh);



  
	// std::cout<<"nlocaldofs="<<decltype(W)::DofsDM::NLocalDofs<<std::endl;
	// auto &dofsdofmap=W_ptr->dofsdofmap();

	// for(Integer i=0;i<dofsdofmap.space_dofs().size();i++)
	// {
	//     std::cout<<"===== space="<<i<<std::endl;
	//     for(Integer j=0;j<dofsdofmap.space_dofs_size(i);j++)
	//     {
	//       // for(Integer k=0;k<sd3[i]->operator[](j).size();k++)
	//          std::cout<<dofsdofmap.space_dofs_get(i,j)<<std::endl;
	//     }

	// }

 // 	for(Integer i=0;i<mesh.points().size();i++)
	// {

	//          std::cout<<mesh.points()[i]<<std::endl;
	// }


	    clock_t begin1 = clock();
		// auto& ordered_entities=entity2dofs.ordered_entities();	     
		clock_t end1 = clock();
	    std::cout<<"TIME BUILD="<< double(end1 - begin1) / CLOCKS_PER_SEC <<std::endl;




 //          

	// for(Integer i=0;i<ordered_entities.size();i++)
	// {
	// 	std::cout<<"level="<<i<<std::endl;
	// 	for(Integer k=0;k<ordered_entities[i].size();k++)
	// 	{
	// 		for(Integer j=0;j<ordered_entities[i][k].size();j++)
	// 		{
	// 			std::cout<< ordered_entities[i][k][j]<< " ";
	// 		}
	// 		std::cout<< std::endl;


	// 	}
	// }
          std::cout<<"END ORDERED ENTITY 2 DOFS"<<std::endl;
          patch_multigrid(solution,A_levels,bL,levels_interp,e2d,3,3,levels.size()-1,100,0.000000001);

         // std::cout<<"START SOLVING PATCH MULTIGRID"<<std::endl;
         // patch_multigrid(solution,A_levels,bL,levels_interp,ordered_entities,5,5,levels.size()-1,50,0.0000000001);

		std::cout<<"END SOLVING PATCH MULTIGRID"<<std::endl;
		std::string output_fileMULTIGRID ="LSFEM_Poisson"+ std::to_string(ManifoldDim) +
		"D_RT" + std::to_string(Order1)+
		"_P" + std::to_string(Order2)+"_outputMULTIGRID.vtk";

		os.close();
		os.open(output_fileMULTIGRID.c_str());
		// write_wtk_isoparametric(os,*context.full_spaces_ptr(),xL,var_namesL,level);
		write_wtk_isoparametric(os,W_ptr,solution,var_namesL,levelL);

	    os.close();

        // std::vector<Integer> a{14,3,2,1,4,2,3,48,3};

        // std::sort(a.begin(),a.end()-1);
        // for(auto it=a.begin();it!=a.end();it++)
        // 	std::cout<<*it<<std::endl;




// 	for(Integer i=0;i<mesh.n_elements();i++)
// 	{
// 	         auto& elem=mesh.elem(i);
// 	         auto& nodes=elem.nodes;
// 	         // std::cout<<"elem == "<< i <<std::endl;
//           //    for(Integer j=0;j<nodes.size();j++)
// 	         // std::cout<<mesh.points()[nodes[j]]<<" ";
//           //    for(Integer j=0;j<nodes.size();j++)
// 	         // std::cout<<nodes[j]<<" ";
// 	         // std::cout<<std::endl;
// 	         if(i==132)
// 	         {
// 	         	std::cout<<"ELEM 132 with level"<<std::endl;
// 	            for(Integer j=0;j<nodes.size();j++)
// 	         std::cout<<nodes[j]<<" ";
         	
// 	         	auto & children=mesh.elem(i).children;
//              for(Integer j=0;j<children.size();j++)
// 	         std::cout<<children[j]<<" ";	         	
// 	         }
// 	}

// auto & children=mesh.elem(132).children;
// 	// for(Integer i=0;i<mesh.points().size();i++)
// 	// {
// 	//          std::cout<<"point == "<< i <<std::endl;
//  //             for(Integer j=0;j<mesh.points()[i].size();j++)
// 	//          std::cout<<mesh.points()[i][j]<<" ";
// 	//          std::cout<<std::endl;
// 	// }


    // auto& bisection=W_ptr->bisection();
 //    auto& tracker=bisection.tracker();


	// std::cout<<"-----dofmap of space="<<0<<std::endl;
	// GetType<typename decltype(W)::DofsDM::ElemDofMap> elemdm;
	// Integer cont1=0;
	// for(Integer i=0;i<dofsdofmap.template dofmap_size<0>();i++)
	// {

	// 		dofsdofmap.template dofmap_get<0>(elemdm,i);
	//          std::cout<<"dofmap == "<< i<<std::endl;
	//   //        std::cout<<elemdm<<std::endl;
	// 	if(!elem_belongs_to_level(mesh,i,level,track)) continue;
			


	// // 	if(i==132)
	// // 	{

	//          auto& elem=mesh.elem(i);
	//          auto& nodes=elem.nodes;
	// 	     dofsdofmap.template dofmap_get<0>(elemdm,i);
	//          std::cout<<"dofmap == "<< i<<std::endl;
	//          std::cout<<elemdm<<std::endl;
	// //          std::cout<<"nodes == "<< i<<std::endl;
 // //             for(Integer k=0;k<nodes.size();k++)
	// //          {
	// //          	std::cout<<"nodes[k] "<<nodes[k]<<std::endl;
	// //          	std::cout<<mesh.points()[nodes[k]]<<" ";
	// //          }


	// // 		for(int j=0;j<children.size();j++)
	// //          {
	// //          auto& elem=mesh.elem(children[j]);
	// //          auto& nodes=elem.nodes;

	// //          dofsdofmap.template dofmap_get<0>(elemdm,children[j]);
	// //          std::cout<<"dofmap == "<< children[j]<<std::endl;
	// //          std::cout<<elemdm<<std::endl;
	// //          std::cout<<"nodes == "<< children[j]<<std::endl;
 // //             for(Integer k=0;k<nodes.size();k++)
	// //          {
	// //          	std::cout<<"nodes[k] "<<nodes[k]<<std::endl;
	// //          	std::cout<<mesh.points()[nodes[k]]<<" ";}

	// //      }
	// //      break;
	// 	// }


	// }

	// std::cout<<"cont1"<<std::endl;
	// std::cout<<cont1<<std::endl;

 
   
 //    auto& signed_normal= mesh.signed_normal().normals();
 //    auto& signs= mesh.signed_normal().sign();



 //    std::cout<<signed_normal[132]<<std::endl;
 //    std::cout<<signed_normal[308]<<std::endl;
 //    std::cout<<signed_normal[293]<<std::endl;
 //    std::cout<<signs[132]<<std::endl;
 //    std::cout<<signs[308]<<std::endl;
 //    std::cout<<signs[293]<<std::endl;



 //    cont1=0;
 //    Integer found=0;
 // 	for(Integer i=0;i<mesh.n_elements();i++)
	// {
	// 	if(!elem_belongs_to_level(mesh,i,level,tracker)) continue;

	// 	if(i==132)
	// 		found=cont1;
	// 	cont1++;

	// }

	// std::cout<<"mesh level "<<level<<" n_elements="<<cont1<<std::endl;

 //    cont1=0;
 //    Integer found1=0;
 //    Integer found2=0;
 // 	for(Integer i=0;i<mesh.n_elements();i++)
	// {
	// 	if(!elem_belongs_to_level(mesh,i,level_finer,tracker)) continue;

	// 	if(i==308)
	// 		found1=cont1;
	// 	if(i==309)
	// 		found2=cont1;
	// 	cont1++;

	// }

	// std::cout<<"mesh level_finer "<<level_finer<<" n_elements="<<cont1<<std::endl;

	// std::cout<<"mesh n_elements = "<< mesh.n_elements()<<std::endl;



	// std::cout<<"mesh level "<<level<<" found="<<found<<std::endl;
	// std::cout<<"mesh level_finer "<<level_finer<<" found1="<<found1<<std::endl;
	// std::cout<<"mesh level_finer "<<level_finer<<" found2="<<found2<<std::endl;


	// std::cout<<"-----dofmap of space="<<1<<std::endl;
	// GetType<typename decltype(W)::DofsDM::ElemDofMap,1> elemdm1;
	// for(Integer i=0;i<dofsdofmap.template dofmap_size<1>();i++)
	// {
	//          dofsdofmap.template dofmap_get<1>(elemdm1,i);
	//          std::cout<<"elem == "<< i<<std::endl;
	//          std::cout<<elemdm1<<std::endl;
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

	// decltype(bilinear_form)::FunctionSpace a1324k(6);
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
