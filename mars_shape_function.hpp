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

namespace mars{





///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////// class SignedNormal: for each face f of the element e                             ///////////////////
////////////// normal[e][f] returns the normal, outward[e][f] returns the sign of the normal    ///////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// returns the sign of val (1,0,-1)
// template <typename T> int sgn(T val) {
//     return (T(0) < val) - (val < T(0));
// }
// returns the order of the indeces of the sorted v
template <typename T>
T argsort(const T &v) {
    // initialize original index locations
    T idx;
    std::iota(idx.begin(), idx.end(), 0);
    // sort indexes based on comparing values in v
    std::sort(idx.begin(), idx.end(),
              [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
    
    return idx;
}
// returns v, resorted based on index 
template <typename T, typename S>
T sort_by_index(const T &v,const S& index) {
    T sortv;
    assert(index.size()==v.size() && "sort_by_index: v and index must have the same length, otherwise last elements in v are not initialized");
    for(Integer ii=0;ii<index.size();ii++)
        sortv[ii]=v[index[ii]];
    return sortv;
}
// returns the sub array of v of indices = index
template <typename S,std::size_t SDim, typename T,std::size_t TDim>
std::array<S,TDim> sub_array(const std::array<S,SDim> &v,const std::array<T,TDim>& index) {
    
    static_assert(TDim<=SDim,"in sub_array the length of the index vector must be smaller than the one of the vector");
    std::array<S,TDim> subvec;
    for(Integer ii=0;ii<TDim;ii++)
        subvec[ii]=v[index[ii]];
    
    return subvec;
}




template<typename Elem>
class SignedNormal;


template<Integer Dim, Integer ManifoldDim>
class SignedNormal<Simplex<Dim,ManifoldDim>>{
private:
    std::vector< Array< Vector< Real, Dim> , ManifoldDim + 1 > >  normal_;
    std::vector< Array<Real, ManifoldDim + 1 > > outward_;
public:
    
    std::vector< Array< Vector< Real, Dim> , ManifoldDim + 1 > > operator () () const { return normal_; };
    // const Vector< Vector< Real, Dim> , ManifoldDim + 1 >& normal(const Integer elem_id) const {return normal_[elem_id];};
    std::vector< Array<Real, ManifoldDim + 1 > > sign() const {return outward_;};
    const Array<Real, ManifoldDim + 1 >& sign(const Integer elem_id) const {return outward_[elem_id];};
    
    template< typename MeshT>
    void init(const MeshT& mesh)
    {
        using Elem=typename MeshT::Elem;
        static_assert(Dim==ManifoldDim, "SignedNormal: computing internal normals of a mesh makes sense only if Dim==ManifoldDim");
        Integer n_elements=mesh.n_elements();
        Integer facenodes_carray[ManifoldDim];
        std::array<Integer,ManifoldDim+1> elemnodes_local;
        std::array<Integer,ManifoldDim+1> elemnodes_global;
        std::array<Integer,ManifoldDim+1> elemnodes_local_sort;
        
        std::array< Integer, ManifoldDim> facenodes_local;
        std::array< Integer, ManifoldDim> facenodes_local_sort;
        std::array<Integer,ManifoldDim> facenodes_global;
        std::array<Integer,ManifoldDim> facenodes_global_sort;
        
        std::array<Vector<Real,Dim>, ManifoldDim+1> points;
        Vector<Vector<Real,Dim>, ManifoldDim> facepoints;
        Vector<Vector<Real,Dim>, ManifoldDim+1> elempoints_sort;
        
        Vector<Real,Dim> facepoint_mean;
        Vector<Real,Dim> elempoint_mean;
        
        // we initialize the simplex
        Simplex<Dim,ManifoldDim> simplex_elem;
        for(Integer nn=0;nn<ManifoldDim+1;nn++)
            simplex_elem.nodes[nn]=nn;
        
        normal_.resize(n_elements);
        outward_.resize(n_elements);
        
        
        // elemnodes_local is the array{ManifoldDim,ManifoldDim-1,...,1,0}
        for(Integer nn=0;nn<ManifoldDim+1 ;nn++)
            elemnodes_local[nn]=nn;
        std::reverse(std::begin(elemnodes_local), std::end(elemnodes_local));
        
        // loop on all the elements
        for(Integer ee=0;ee<mesh.n_elements();ee++)
        {
            Elem elem=mesh.elem(ee);
            elemnodes_global=elem.nodes;
            
            for(Integer mm=0;mm<points.size();mm++)
                points[mm]=mesh.point(elemnodes_global[mm]);
            
            // loop on all the faces of the simplex
            // for Simplex={0,1,2,3}, the order is: 0-1-2, 0-1-3, 0-2-3, 1-2-3
            for(Integer mm=0;mm<ManifoldDim+1;mm++)
            {
                // facenodes_local is a std::array containing the local nodes of the face mm
                Combinations<ManifoldDim + 1,ManifoldDim>::generate(mm,facenodes_carray);
                std::copy(std::begin(facenodes_carray), std::end(facenodes_carray), std::begin(facenodes_local));
                // we reorder the local nodes based on the sorting order of the corresponding global nodes
                facenodes_global=sub_array(elemnodes_global,facenodes_local);
                facenodes_global_sort=argsort(facenodes_global);
                facenodes_local_sort=sort_by_index(facenodes_local,facenodes_global_sort);
                
                // in elemnodes_local_sort we have facenodes_local_sort and the remaining node, elemnodes_local[mm]
                for(Integer nn=0;nn<ManifoldDim;nn++)
                    elemnodes_local_sort[nn]=facenodes_local_sort[nn];
                
                elemnodes_local_sort[ManifoldDim]=elemnodes_local[mm];
                
                for(Integer nn=0;nn<ManifoldDim+1;nn++)
                    elempoints_sort[nn]=points[elemnodes_local_sort[nn]];
                
                // we create the face midpoint, the element midpoint and their difference
                for(Integer nn=0;nn<ManifoldDim;nn++)
                    facepoints[nn]=points[facenodes_local[nn]];
                
                facepoint_mean = facepoints.Tmean();
                elempoint_mean=elempoints_sort.Tmean();
                auto diff=facepoint_mean-elempoint_mean;
                
                auto n = normal(simplex_elem, elempoints_sort);
                outward_[ee][mm]=Sign(dot(diff,n));
                // save the normal and if its outward or inward
                // if the face is on the boundary, then it must be outward
                if(elem.side_tags[mm]>=0 && outward_[ee][mm]==-1)
                {
                    // if the normal ins inward, force it to be outward
                    normal_[ee][mm]=-n;
                    outward_[ee][mm]=1;
                }
                else
                {
                    normal_[ee][mm]=n;
                }
            }
        }
    };
    
    SignedNormal(){}
    
    template< typename MeshT>
    SignedNormal(const MeshT& mesh)
    {init(mesh);};
    
    template< typename MeshT>
    void print(const MeshT& mesh)
    {
        using Elem=typename MeshT::Elem;
        Vector<Vector<Real,Dim>, ManifoldDim+1> elempoints_sort;
        std::array<Integer,ManifoldDim+1> nodes;
        std::vector<Vector<Real,Dim>> points;
        Integer facenodes_carray[ManifoldDim];
        std::array< Integer, ManifoldDim> facenodes_local;
        
        for(Integer ee=0;ee<normal_.size();ee++)
        {
            std::cout<<std::endl<<"ELEMENT ID == "<<ee<<" WITH NODES"<<std::endl;
            Elem elem=mesh.elem(ee);
            nodes=elem.nodes;
            for(Integer mm=0;mm<ManifoldDim+1;mm++)
                std::cout<< nodes[mm] <<" ";
            std::cout<<std::endl;
            
            for(Integer mm=0;mm<ManifoldDim+1;mm++)
            {
                Combinations<ManifoldDim + 1,ManifoldDim>::generate(mm,facenodes_carray);
                std::copy(std::begin(facenodes_carray), std::end(facenodes_carray), std::begin(facenodes_local));
                std::cout<<"facenodes_local== "<<std::endl;
                for(Integer ii=0;ii<ManifoldDim;ii++)
                    std::cout<<facenodes_local[ii]<<" ";
                std::cout<<"normal== ";
                normal_[ee][mm].describe(std::cout);
                std::cout<<"outward/inward(+1/-1)== "<<outward_[ee][mm]<<std::endl;
            }
        }
    };
    
    
};











































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

 static constexpr auto dofs_aux=std::tuple_cat(Args::dofs()...);

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

template<typename Elem,typename Operator,Integer FEFamily,Integer Order>
class DofsPoints;


template<typename DofsPoints_>
class DofsPointsType;

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












template<Integer Dim,Integer ManifoldDim>
class ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, IdentityOperator, LagrangeFE, 0>
{
 public: 
  using Output=Vector<Matrix<Real, 1, 1>,1>;
constexpr inline static void 
apply(const Vector<Real,ManifoldDim>& point, Output & func)
{
    Output func2(1);       
    func=func2;
}
};


template<Integer Dim>
class ReferenceShapeFunctionValue<Simplex<Dim,1>, IdentityOperator, LagrangeFE, 1>
{
 public: 
  using Output=Vector<Matrix<Real, 1, 1>,2>;
constexpr inline static void 
apply(const Vector<Real,1>& point, Output & func)
{
    Output func2((1. - point[0]), // 1 in (0)
                                       point[0]);       // 1 in (1)
    func=func2;
}
};

template<Integer Dim>
class ReferenceShapeFunctionValue<Simplex<Dim,2>, IdentityOperator, LagrangeFE, 1>
{
 public: 
  using Output=Vector<Matrix<Real, 1, 1>,3>;
constexpr inline static void 
apply(const Vector<Real,2>& point, Output & func)
{
    Output func2((1. - point[0] - point[1]), // 1 in (0,0)
                  point[0],                  // 1 in (1,0)
                  point[1]);                 // 1 in (0,1)
    func=func2;
}
};


// template<typename Elem,typename Operator,Integer FEFamily,Integer Order,typename Output,typename Point>
// constexpr void value(const Point& point,Output& output);

// template<Integer Dim>
// constexpr void value<Simplex<Dim,2>, IdentityOperator, LagrangeFE, 1>
//  (const Vector<Real,2>& point, Vector<Matrix<Real, 1, 1>,3> & func)
// {
//     Vector<Matrix<Real, 1, 1>,3> func2((1. - point[0] - point[1]), // 1 in (0,0)
//                                        point[0],                  // 1 in (1,0)
//                                        point[1]);                 // 1 in (0,1)
//     func=func2;
// }



template<Integer Dim>
class ReferenceShapeFunctionValue<Simplex<Dim,2>, GradientOperator, LagrangeFE, 1>
{
 public: 
  using Output=Vector<Matrix<Real, 2, 1>,3>;
 constexpr inline static void 
 apply(const Vector<Real,2>& point, Output & func)
  {
      const auto& xi=point[0];
      const auto& eta=point[1];
      const Real zeta = 1. - xi - eta;
      Output func2({-1,-1},
                                         // Vector<Vector<Real, 2>,3> func2({-1,-1},
                                         {+1, 0},
                                         { 0,+1});
      func=func2;
  }
};

// template<Integer Dim>
// constexpr void value<Simplex<Dim,2>, GradientOperator, LagrangeFE, 1>
// (const Vector<Real,2>& point, Vector<Matrix<Real, 2, 1>,3> & func)
// {
//     const auto& xi=point[0];
//     const auto& eta=point[1];
//     const Real zeta = 1. - xi - eta;
//     Vector<Matrix<Real, 2, 1>,3> func2({-1,-1},
//                                        // Vector<Vector<Real, 2>,3> func2({-1,-1},
//                                        {+1, 0},
//                                        { 0,+1});
//     func=func2;
// }
template<Integer Dim>
class ReferenceShapeFunctionValue<Simplex<Dim,2>, TraceOperator, LagrangeFE, 1>
{
 public: 
  using Output=Vector<Matrix<Real, 1, 1>,2>;
constexpr inline static void 
apply(const Vector<Real,1>& point, Output & func)
  {
      ReferenceShapeFunctionValue<Simplex<Dim,1>, IdentityOperator, LagrangeFE, 1>::apply(point,func);
  }
};








template<Integer Dim>
class ReferenceShapeFunctionValue<Simplex<Dim,2>, IdentityOperator, LagrangeFE, 2>
{
 public:
  using Output=Vector<Matrix<Real, 1, 1>,6>;
constexpr inline static void 
apply(const Vector<Real,2>& point, Output & func)
{
    const auto& xi=point[0];
    const auto& eta=point[1];
    const Real zeta = 1. - xi - eta;
    Output func2(2.*zeta*(zeta-0.5), // 1 in (0,0)
                                       2.*xi*(xi-0.5),     // 1 in (1,0)
                                       2.*eta*(eta-0.5),   // 1 in (0,1)
                                       4.*zeta*xi,         // 1 in (0.5,0)
                                       4.*eta*zeta,        // 1 in (0,0.5)
                                       4.*xi*eta);         // 1 in (0.5,0.5)
    func=func2;
}
};


// template<Integer Dim>
// constexpr void value<Simplex<Dim,2>, IdentityOperator, LagrangeFE, 2>
//  (const Vector<Real,2>& point, Vector<Matrix<Real, 1, 1>,6> & func)
// {
//     const auto& xi=point[0];
//     const auto& eta=point[1];
//     const Real zeta = 1. - xi - eta;
//     Vector<Matrix<Real, 1, 1>,6> func2(2.*zeta*(zeta-0.5), // 1 in (0,0)
//                                        2.*xi*(xi-0.5),     // 1 in (1,0)
//                                        2.*eta*(eta-0.5),   // 1 in (0,1)
//                                        4.*zeta*xi,         // 1 in (0.5,0)
//                                        4.*eta*zeta,        // 1 in (0,0.5)
//                                        4.*xi*eta);         // 1 in (0.5,0.5)
//     func=func2;
// }

template<Integer Dim>
class ReferenceShapeFunctionValue<Simplex<Dim,2>, GradientOperator, LagrangeFE, 2>
{
public:
  using Output=Vector<Matrix<Real, 2, 1>,6>;
constexpr inline static void 
apply(const Vector<Real,2>& point, Output & func)
{
    const auto& xi=point[0];
    const auto& eta=point[1];
    const Real zeta = 1. - xi - eta;
    const Real dxi_dxi    = 1.;
    const Real deta_deta  = 1.;
    const Real dzeta_dxi  = -1.;
    const Real dzeta_deta = -1.;
    // Vector<Vector<Real, 2>,6> func2(
    Output func2(
                                       {2.*zeta*dzeta_dxi  + 2*dzeta_dxi *(zeta-0.5), 2.*zeta*dzeta_deta + 2*dzeta_deta*(zeta-0.5)},
                                       {2.*xi*dxi_dxi  + 2.*dxi_dxi *(xi-0.5),0},
                                       {0,2.*eta*deta_deta + 2.*deta_deta*(eta-0.5)},
                                       {4.*zeta*dxi_dxi+4.*dzeta_dxi*xi,4.*dzeta_deta*xi},
                                       {4.*dxi_dxi*eta,4.*xi*deta_deta},
                                       {4.*eta*dzeta_dxi,4.*eta*dzeta_deta + 4.*deta_deta*zeta});
    
    func=func2;
  }
};





// template<Integer Dim>
// constexpr void value<Simplex<Dim,2>, GradientOperator, LagrangeFE, 2>
// // (const Vector<Real,2>& point, Vector<Vector<Real, 2>,6> & func)
//  (const Vector<Real,2>& point, Vector<Matrix<Real, 2, 1>,6> & func)
// {
//     const auto& xi=point[0];
//     const auto& eta=point[1];
//     const Real zeta = 1. - xi - eta;
//     const Real dxi_dxi    = 1.;
//     const Real deta_deta  = 1.;
//     const Real dzeta_dxi  = -1.;
//     const Real dzeta_deta = -1.;
//     // Vector<Vector<Real, 2>,6> func2(
//     Vector<Matrix<Real, 2, 1>,6> func2(
//                                        {2.*zeta*dzeta_dxi  + 2*dzeta_dxi *(zeta-0.5), 2.*zeta*dzeta_deta + 2*dzeta_deta*(zeta-0.5)},
//                                        {2.*xi*dxi_dxi  + 2.*dxi_dxi *(xi-0.5),0},
//                                        {0,2.*eta*deta_deta + 2.*deta_deta*(eta-0.5)},
//                                        {4.*zeta*dxi_dxi+4.*dzeta_dxi*xi,4.*dzeta_deta*xi},
//                                        {4.*dxi_dxi*eta,4.*xi*deta_deta},
//                                        {4.*eta*dzeta_dxi,4.*eta*dzeta_deta + 4.*deta_deta*zeta});
    
//     func=func2;
// }



template<Integer Dim>
class ReferenceShapeFunctionValue<Simplex<Dim,2>, IdentityOperator, RaviartThomasFE, 0>
{
public:
  using Output=Vector<Matrix<Real, 2, 1>,3>;
 constexpr inline static void 
 apply(const Vector<Real,2>& point, Output & func)
{
    const auto& xi=point[0];
    const auto& eta=point[1];
    Output func2{{xi,eta-1},{xi-1,eta},{xi,eta}};
    func=func2;
}

};

//shape_function_coefficients_init for: Simplex<2,2>, LagrangeFE, 2
// do nothing: lagrange shape functions do not need any coefficient

// template<Integer Dim>
// constexpr void value<Simplex<Dim,2>, IdentityOperator, RaviartThomasFE, 0>
//  (const Vector<Real,2>& point, Vector<Matrix<Real, 2, 1>,3> & func)
// {
//     const auto& xi=point[0];
//     const auto& eta=point[1];
//     Vector<Matrix<Real, 2, 1>,3> func2{{xi,eta-1},{xi-1,eta},{xi,eta}};
//     func=func2;
// }

template<Integer Dim>
class ReferenceShapeFunctionValue<Simplex<Dim,2>, DivergenceOperator, RaviartThomasFE, 0>
{
public:
  using Output=Vector<Matrix<Real, 1, 1>,3>;
constexpr inline static void 
apply(const Vector<Real,2>& point, Output & func)
{
    Output func2{{2},{2},{2}};
    func=func2;
}
};



// template<Integer Dim>
// constexpr void value<Simplex<Dim,2>, DivergenceOperator, RaviartThomasFE, 0>
//  (const Vector<Real,2>& point, Vector<Matrix<Real, 1, 1>,3> & func)
// {
//     Vector<Matrix<Real, 1, 1>,3> func2{{2},{2},{2}};
//     func=func2;
// }

template<Integer Dim>
class ReferenceShapeFunctionValue<Simplex<Dim,2>, IdentityOperator, RaviartThomasFE, 1>
{
public:
  using Output=Vector<Matrix<Real, 2, 1>,8>;
constexpr inline static void 
apply (const Vector<Real,2>& point, Output & func)
{
    const auto& xi=point[0];
    const auto& eta=point[1];
    
    Output func2{
        {(1. - xi - eta)*xi,(1. - xi - eta)*(eta-1)},   // 0 in (1,0), (0,1), non-zero normal on edge0
        {xi*xi,xi*(eta-1)},                             // 0 in (0,0), (0,1), non-zero normal on edge0
        {(1. - xi - eta)*(xi-1),(1. - xi - eta)*(eta)}, // 0 in (1,0), (0,1), non-zero normal on edge1
        {eta*(xi-1),eta*eta},                           // 0 in (0,0), (1,0), non-zero normal on edge1
        {xi*xi,xi*eta},                                 // 0 in (0,0), (0,1), non-zero normal on edge2
        {eta*xi,eta*eta},                               // 0 in (0,0), (1,0), non-zero normal on edge2
        {eta*xi,eta*(eta-1)},                           // normal 0 on all edges, element-dof
        {xi*(xi-1),xi*eta}                              // normal 0 on all edges, element-dof
    };
    func=func2;
}
};
// template<Integer Dim>
// constexpr void value<Simplex<Dim,2>, IdentityOperator, RaviartThomasFE, 1>
//  (const Vector<Real,2>& point, Vector<Matrix<Real, 2, 1>,8> & func)
// {
//     const auto& xi=point[0];
//     const auto& eta=point[1];
    
//     Vector<Matrix<Real, 2, 1>,8> func2{
//         {(1. - xi - eta)*xi,(1. - xi - eta)*(eta-1)},   // 0 in (1,0), (0,1), non-zero normal on edge0
//         {xi*xi,xi*(eta-1)},                             // 0 in (0,0), (0,1), non-zero normal on edge0
//         {(1. - xi - eta)*(xi-1),(1. - xi - eta)*(eta)}, // 0 in (1,0), (0,1), non-zero normal on edge1
//         {eta*(xi-1),eta*eta},                           // 0 in (0,0), (1,0), non-zero normal on edge1
//         {xi*xi,xi*eta},                                 // 0 in (0,0), (0,1), non-zero normal on edge2
//         {eta*xi,eta*eta},                               // 0 in (0,0), (1,0), non-zero normal on edge2
//         {eta*xi,eta*(eta-1)},                           // normal 0 on all edges, element-dof
//         {xi*(xi-1),xi*eta}                              // normal 0 on all edges, element-dof
//     };
//     func=func2;
// }

template<Integer Dim>
class ReferenceShapeFunctionValue<Simplex<Dim,2>, DivergenceOperator, RaviartThomasFE, 1>
{
public:
  using Output=Vector<Matrix<Real, 1, 1>,8>;
constexpr inline static void 
apply(const Vector<Real,2>& point, Output & func)
{
    const auto& xi=point[0];
    const auto& eta=point[1];
    Output func2{
        {3*(1-xi-eta)},
        {3*xi},
        {3*(1-xi-eta)},
        {3*eta},
        {3*xi},
        {3*eta},
        {3*eta},
        {3*xi}
    };
    func=func2;
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
      ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim-1>, IdentityOperator, LagrangeFE, Order>::apply(point,func);
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











template<Integer Dim>
class SingleShapeFunctionCoefficientsCollection<Simplex<Dim,2>, RaviartThomasFE, 0>
{
 public: 
constexpr inline static  void 
apply(const Array<Real, 3 >& outward,Array<Real, 3 >& coeff)
{
    coeff[0]=outward[0];
    coeff[1]=outward[1];
    coeff[2]=outward[2];
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



template<>
void shape_function_coefficients_init<Simplex<2,2>, RaviartThomasFE, 0>
 (const Vector<Real, 3 >& outward,Vector<Real, 3 >& coeff)
{
    coeff[0]=outward[0];
    coeff[1]=outward[1];
    coeff[2]=outward[2];
}





template<>
void shape_function_coefficients_init<Simplex<2,2>, RaviartThomasFE, 1>
 (const Vector<Real, 3 >& outward,Vector<Real, 8 >& coeff)
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
        // decltype(func_values_) ok1(1);
        // SingleType ok2(2);
        // std::cout<<"NComponents="<<NComponents<<std::endl;
        // std::cout<<"NQPoints="<<NQPoints<<std::endl;
        
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
                }
                
            }
        }
        
        std::cout<<"ShapeFunction non trace end"<<std::endl;
    }
    
    void init(const Array<Real,Ndofs> &alpha)
    {
        std::cout<<"init RT elements (coeffs)"<<std::endl;
        
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
                    std::cout<< "n_dof, qp, n_tot_, n_comp, n_=("<<n_dof<<", "<<qp<<", "<< n_tot_<<", "<<n_comp<<", "<< n_<<")"<<std::endl;
                    std::cout<< "func_values_="<<func_values_[n_tot_][qp]<<std::endl;
                    func_tmp_=alpha[n_dof] * mapping * reference_values[n_dof][qp];
                    std::cout<< "func_tmp_="<<func_tmp_<<std::endl;
                    assign<NComponents>(func_values_[n_tot_][qp],func_tmp_,n_,0);
                    std::cout<< "func_values_ after="<<func_values_[n_tot_][qp]<<std::endl;
                }
                
            }
        }
        std::cout<<"init end"<<std::endl;
        
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
    
    using trace_type=typename decltype(trace)::value_type;
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
        
     std::cout<<"ShapeFunction  trace end"<<std::endl;

    }
    
    void init(const Array<Real,NdofsVolume> &beta, const Integer face)
    {

        //face=0
        // const Integer face=0;
        // auto alpha=subarray(beta,trace[face]);
        subarray(alpha_,beta,trace[face]);
        
        std::cout<<"init RT elements (coeffs)"<<std::endl;
        
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
                    std::cout<< "n_dof, qp, n_tot_, n_comp, n_=("<<n_dof<<", "<<qp<<", "<< n_tot_<<", "<<n_comp<<", "<< n_<<")"<<std::endl;
                    std::cout<< "func_values_="<<func_values_[n_tot_][qp]<<std::endl;
                    func_tmp_=alpha_[n_dof] * mapping * reference_values[n_dof][qp];
                    std::cout<< "func_tmp_="<<func_tmp_<<std::endl;
                    assign<NComponents>(func_values_[n_tot_][qp],func_tmp_,n_,0);
                    std::cout<< "func_values_ after="<<func_values_[n_tot_][qp]<<std::endl;
                }
                
            }
        }
        std::cout<<"init end"<<std::endl;
        
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






