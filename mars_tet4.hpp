#ifndef MARS_TET4
#define MARS_TET4

#include "mars_gauss_points.hpp"
#include "mars_affine_transformation.hpp"
#include "generation/mars_mesh_kokkos.hpp"
namespace mars {

  template<typename T, int N_>
  class Array {
  public:
    constexpr static const int N = N_; 
    MARS_INLINE_FUNCTION constexpr Array()
    {}

    MARS_INLINE_FUNCTION constexpr int size() const
    {
        return N;
    }

    T values[N];
  };

  class Tet4 {
  public:

    MARS_INLINE_FUNCTION constexpr Tet4()
    {}


    MARS_INLINE_FUNCTION constexpr int size() const
    {
        return 4;
    }

   /*evaluation basis function*/
    MARS_INLINE_FUNCTION static void eval_phi(const Real *p, Array<Real, 4> &ret) {
        ret.values[0] = 1.0 - p[0] - p[1] - p[2];
        ret.values[1] = p[0];
        ret.values[2] = p[1];
        ret.values[3] = p[2];
    }


    /*evaluation gradient of basis function*/
    MARS_INLINE_FUNCTION static void eval_grad(const Real *p, Array<Real, 12> &ret) {
        
      //fn 0
        ret.values[0] = -1.0;
        ret.values[1] = -1.0;
        ret.values[2] = -1.0;
      //fn 1
        ret.values[3] =  1.0;
        ret.values[4] =  0.0;
        ret.values[5] =  0.0;
      //fn 2
        ret.values[6] =  0.0;
        ret.values[7] =  1.0;
        ret.values[8] =  0.0;
      //fn 3
        ret.values[9]  = 0.0;
        ret.values[10] = 0.0;
        ret.values[11] = 1.0;
    }

   MARS_INLINE_FUNCTION Real static measure ()
   {
     return 1.0/6;
   }

private:
    // MARS_INLINE_FUNCTION constexpr Array<Array<Real, 4>, QGauss<3, 2>::N> init(const QGauss<3, 2> &q)
    // {
    //     static_assert(QGauss<3, 2>::N == 4, "expected quad rule with 4 values");
    //     Array<Array<Real, 4>, QGauss<3, 2>::N> values;
    //     values.values[0] = eval(q.point(0));
    //     values.values[1] = eval(q.point(1));
    //     values.values[2] = eval(q.point(2));
    //     values.values[3] = eval(q.point(3));
    //     return values;
    // }


  };

  template<class Q, class FE>
  class RefMassMatrix {};

  template<>
  class RefMassMatrix<QGauss<3, 2>, Tet4> {
  public:

    static constexpr const Real c = 0.1;

    MARS_INLINE_FUNCTION constexpr RefMassMatrix()
    : value{c}

    {}

    MARS_INLINE_FUNCTION constexpr const Real elem_value() const
    {
        return value;
    }


     Real value;

  };

  template<class Q, class FE>
  class AssembleMassMatrix {
  public:
        ParallelMesh3 pMesh;
        Kokkos::View<double*[4][4]> element_matrix;

        AssembleMassMatrix(
            ParallelMesh3 pMesh,
            Kokkos::View<double*[4][4]> element_matrix)
        : pMesh(pMesh), element_matrix(element_matrix) {}

        AssembleMassMatrix() {}



        MARS_INLINE_FUNCTION void operator()(int index) const
        {
           
              auto e = pMesh.elem(index);

              constexpr const Q q;

              constexpr const FE phi;

              RefMassMatrix<Q, FE> a;

              AffineTransform<3,3> Map;

              Real Jacobian[3*3];

              Real J; 





              Map.affine_jacobian(pMesh, e, Jacobian, J);

              std::cout<<"value j"<<J<<std::endl;

              for(int i=0; i<4; i++)
                  for(int j=0; j<4; j++)
                    element_matrix(index,i,j) = a.elem_value() * J * phi.measure();

              std::cout<<"outside op"<<std::endl;

            
        }


     
      


            // Kokkos::parallel_for(q.N, KOKKOS_LAMBDA ( int k ){

            // dx[k] = detJ * phi.measure();

            // });         
 

 
    };


  /*template<typename Scalar>
  KERNEL_PREFIX void gradients(const Scalar* x, Scalar* values_per_fn)
  {

  //fn 0
    values_per_fn[0] = -1.0;
    values_per_fn[1] = -1.0;
    values_per_fn[2] = -1.0;
  //fn 1
    values_per_fn[3] =  1.0;
    values_per_fn[4] =  0.0;
    values_per_fn[5] =  0.0;
  //fn 2
    values_per_fn[6] =  0.0;
    values_per_fn[7] =  1.0;
    values_per_fn[8] =  0.0;
  //fn 3
    values_per_fn[9]  = 0.0;
    values_per_fn[10] = 0.0;
    values_per_fn[11] = 1.0;
  }


  template<typename Scalar>
  KERNEL_PREFIX Scalar determinant3x3(const Scalar* J)
  {
  //hardwired "3x3" in function-name allows us to assume that J has length 9:

    Scalar J00 = J[0];
    Scalar J01 = J[1];
    Scalar J02 = J[2];

    Scalar J10 = J[3];
    Scalar J11 = J[4];
    Scalar J12 = J[5];

    Scalar J20 = J[6];
    Scalar J21 = J[7];
    Scalar J22 = J[8];

    Scalar term0 = J22*J11 - J21*J12;
    Scalar term1 = J22*J01 - J21*J02;
    Scalar term2 = J12*J01 - J11*J02;

    Scalar detJ = J00*term0 - J10*term1 + J20*term2;

    return detJ;
  }


  template<typename Scalar>
  KERNEL_PREFIX void gradients_and_detJ(const Scalar* elemNodeCoords,
                                        const Scalar* grad_vals,
                                        Scalar& detJ)
  {

    const Scalar zero = 0;

    Scalar J00 = zero;
    Scalar J01 = zero;
    Scalar J02 = zero;

    Scalar J10 = zero;
    Scalar J11 = zero;
    Scalar J12 = zero;

    Scalar J20 = zero;
    Scalar J21 = zero;
    Scalar J22 = zero;

    size_t i_X_spatialDim = 0;
    for(size_t i=0; i<numNodesPerElem; ++i) {
      J00 += grad_vals[i_X_spatialDim+0]*elemNodeCoords[i_X_spatialDim+0];
      J01 += grad_vals[i_X_spatialDim+0]*elemNodeCoords[i_X_spatialDim+1];
      J02 += grad_vals[i_X_spatialDim+0]*elemNodeCoords[i_X_spatialDim+2];

      J10 += grad_vals[i_X_spatialDim+1]*elemNodeCoords[i_X_spatialDim+0];
      J11 += grad_vals[i_X_spatialDim+1]*elemNodeCoords[i_X_spatialDim+1];
      J12 += grad_vals[i_X_spatialDim+1]*elemNodeCoords[i_X_spatialDim+2];

      J20 += grad_vals[i_X_spatialDim+2]*elemNodeCoords[i_X_spatialDim+0];
      J21 += grad_vals[i_X_spatialDim+2]*elemNodeCoords[i_X_spatialDim+1];
      J22 += grad_vals[i_X_spatialDim+2]*elemNodeCoords[i_X_spatialDim+2];

      i_X_spatialDim += spatialDim;
    }

    Scalar term0 = J22*J11 - J21*J12;
    Scalar term1 = J22*J01 - J21*J02;
    Scalar term2 = J12*J01 - J11*J02;

    detJ = J00*term0 - J10*term1 + J20*term2;
 }

  template<typename Scalar>
  KERNEL_PREFIX void gradients_and_invJ_and_detJ(const Scalar* elemNodeCoords,
                                                 const Scalar* grad_vals,
                                                 Scalar* invJ,
                                                 Scalar& detJ)
  {

    const Scalar zero = 0;

    Scalar J00 = zero;
    Scalar J01 = zero;
    Scalar J02 = zero;

    Scalar J10 = zero;
    Scalar J11 = zero;
    Scalar J12 = zero;

    Scalar J20 = zero;
    Scalar J21 = zero;
    Scalar J22 = zero;

    size_t i_X_spatialDim = 0;
    for(size_t i=0; i<numNodesPerElem; ++i) {
      J00 += grad_vals[i_X_spatialDim+0]*elemNodeCoords[i_X_spatialDim+0];
      J01 += grad_vals[i_X_spatialDim+0]*elemNodeCoords[i_X_spatialDim+1];
      J02 += grad_vals[i_X_spatialDim+0]*elemNodeCoords[i_X_spatialDim+2];

      J10 += grad_vals[i_X_spatialDim+1]*elemNodeCoords[i_X_spatialDim+0];
      J11 += grad_vals[i_X_spatialDim+1]*elemNodeCoords[i_X_spatialDim+1];
      J12 += grad_vals[i_X_spatialDim+1]*elemNodeCoords[i_X_spatialDim+2];

      J20 += grad_vals[i_X_spatialDim+2]*elemNodeCoords[i_X_spatialDim+0];
      J21 += grad_vals[i_X_spatialDim+2]*elemNodeCoords[i_X_spatialDim+1];
      J22 += grad_vals[i_X_spatialDim+2]*elemNodeCoords[i_X_spatialDim+2];

      i_X_spatialDim += spatialDim;
    }

    Scalar term0 = J22*J11 - J21*J12;
    Scalar term1 = J22*J01 - J21*J02;
    Scalar term2 = J12*J01 - J11*J02;

    detJ = J00*term0 - J10*term1 + J20*term2;

    Scalar inv_detJ = 1.0/detJ;

    invJ[0] =  term0*inv_detJ;
    invJ[1] = -term1*inv_detJ;
    invJ[2] =  term2*inv_detJ;

    invJ[3] = -(J22*J10 - J20*J12)*inv_detJ;
    invJ[4] =  (J22*J00 - J20*J02)*inv_detJ;
    invJ[5] = -(J12*J00 - J10*J02)*inv_detJ;

    invJ[6] =  (J21*J10 - J20*J11)*inv_detJ;
    invJ[7] = -(J21*J00 - J20*J01)*inv_detJ;
    invJ[8] =  (J11*J00 - J10*J01)*inv_detJ;
   }*/



          // for(std::size_t k = 0; k < n_qp_; ++k) {
          //     dx_[k] = det_JxM * q.weights[k];
          // }

         

   
            //read the element points
            //get J
            //
            // Real vals[4 * 4];
            
            // Array<Real, 4> ret; 

            // phi.eval(q.points, ret);

            // for(int i=0; i<phi.size(); i++)
            // {
                
            //   for(int j=0; j<phi.size(); j++)
            //   {

            //     for(int k=0; k<q.N; k++)
            //     {               
            //         std::cout<<"ret.values["<<i<<"]"<<"["<<k<<"]"<<ret.values[i,k]<<std::endl;

            //         vals[i, j] += q.weights[k] * 1.0 * ret.values[i,k] * ret.values[j,k];
            //     }

            //   }
            // }

            // std::cout<<vals[0,0]<<std::endl;
            // std::cout<<vals[0,1]<<std::endl;
            // std::cout<<vals[0,2]<<std::endl;
            // std::cout<<vals[0,3]<<std::endl;
            // std::cout<<vals[1,0]<<std::endl;
            // std::cout<<vals[1,1]<<std::endl;
            // std::cout<<vals[1,2]<<std::endl;
            // std::cout<<vals[1,3]<<std::endl;
            // std::cout<<vals[2,0]<<std::endl;
            // std::cout<<vals[2,1]<<std::endl;
            // std::cout<<vals[2,2]<<std::endl;
            // std::cout<<vals[2,3]<<std::endl;
            // std::cout<<vals[3,0]<<std::endl;
            // std::cout<<vals[3,1]<<std::endl;
            // std::cout<<vals[3,2]<<std::endl;
            // std::cout<<vals[3,3]<<std::endl;

}
#endif
