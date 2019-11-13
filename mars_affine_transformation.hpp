#ifndef MARS_AFFINETRANSFORM
#define MARS_AFFINETRANSFORM
#include "generation/mars_static_math_kokkos.hpp"

#include <map>
#include <vector>
#include <memory>
#include <cassert>
#include <set>

#include "generation/mars_edge_kokkos.hpp"
#include "generation/mars_edge_select_kokkos.hpp"
#include "mars_longest_edge.hpp"
#include "mars_fwd.hpp"
#include "generation/mars_mesh_kokkos.hpp"

namespace mars {
template<Integer Dim, Integer ManifoldDim>
class AffineTransform {
  public:

    MARS_INLINE_FUNCTION constexpr AffineTransform()
    {}



    /*for gpu and cpu*/
    MARS_INLINE_FUNCTION static void affine_jacobian(ParallelMesh3 pMesh, const Simplex<Dim, ManifoldDim, KokkosImplementation> &e, Real J[3*3], Real& detJ) {
          
          
          const auto &v0 = e.nodes[0];

          auto p0 = pMesh.point(v0);

          std::cout<<"INSIDE"<<e.n_nodes()<<std::endl;

          std::cout<<"INSIDE p0.x=>"<<p0[0]<<std::endl;
          std::cout<<"INSIDE p0.y=>"<<p0[1]<<std::endl;
          std::cout<<"INSIDE p0.z=>"<<p0[2]<<std::endl;

          for (int j=0; j<Dim; j++ )
          {
            
            const int i_offset = j * Dim;

            for ( int i = 1; i < Dim+1; ++i ) 
            {
               const auto &v = e.nodes[i];

               const auto p = pMesh.point(v);
               std::cout<<"INSIDE p0.x=>"<<p[0]<<std::endl;
               std::cout<<"INSIDE p0.y=>"<<p[1]<<std::endl;
               std::cout<<"INSIDE p0.z=>"<<p[2]<<std::endl; 
                   /*mat[0]=p1x - p0x, mat[1]=p2x - p0x, mat[3]=p3x - p0x*/

               J[i_offset + i-1]= p[j]-p0[j];
               std::cout<<"J"<< J[i_offset + i-1]<<std::endl;
             }
          }

            Real J00=J[0];
            Real J01=J[1];
            Real J02=J[2];

            Real J10=J[3];
            Real J11=J[4];
            Real J12=J[5];

            Real J20=J[6];
            Real J21=J[7];
            Real J22=J[8];


            detJ = J[0] * J[4] * J[8]  +
                   J[1] * J[5] * J[6]  +
                   J[2] * J[3] * J[7]  -
                   J[0] * J[5] * J[7]  -
                   J[1] * J[3] * J[8]  -
                   J[2] * J[4] * J[6];

        std::cout<<"OUTSIDE"<<v0<<std::endl;
               
      }
  };
}
#endif 


    // template<Integer Dim, Integer ManifoldDim>
    // inline void jacobian(const Simplex<Dim, ManifoldDim>  &simplex,
    //                      const std::vector<Vector<Real, Dim>> &points,
    //                      Matrix<Real, Dim, ManifoldDim> &J)
    // {
    //     static_assert(Dim >= ManifoldDim, "Dim must be greater or equal ManifoldDim");
        
    //     auto n = n_nodes(simplex);
        
    //     Vector<Real, Dim> v0 = points[simplex.nodes[0]];
        
    //     for(Integer i = 1; i < n; ++i) {
    //         const auto &vi = points[simplex.nodes[i]];
    //         J.col(i-1, vi - v0);
    //     }
    // }