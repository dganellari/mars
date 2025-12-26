#ifndef FE_HEX_H
#define FE_HEX_H


// Only include Eigen on host, not in CUDA device code
#ifndef __CUDA_ARCH__
#include <Eigen/Dense>
#include <Eigen/Core>
#endif

template <typename scalar, typename ExecSpace = Kokkos::DefaultExecutionSpace>
struct FEhex
{
    static constexpr unsigned dim = 3;
    static constexpr unsigned nodesPerElement = 8;
    static constexpr unsigned numScsIp = 12;
    static constexpr unsigned numScvIp = 8;

    // views for shape function data
    using HexDerivView = Kokkos::View<scalar[numScsIp][nodesPerElement][dim],
                                      Kokkos::LayoutLeft,
                                      ExecSpace>;
    using HexDerivViewHost = typename HexDerivView::HostMirror;
    using ConstHexDerivView =
        Kokkos::View<const scalar[numScsIp][nodesPerElement][dim],
                     Kokkos::LayoutLeft,
                     ExecSpace>;

    using SFView = Kokkos::
        View<scalar[numScsIp][nodesPerElement], Kokkos::LayoutLeft, ExecSpace>;
    using SFViewHost = typename SFView::HostMirror;
    using ConstSFView = Kokkos::View<const scalar[numScsIp][nodesPerElement],
                                     Kokkos::LayoutLeft,
                                     ExecSpace>;
    // standard integration location
    static constexpr double intgLoc[numScsIp * dim] = {
        0.00,  -0.25, -0.25, // surf 1    1->2
        0.25,  0.00,  -0.25, // surf 2    2->3
        0.00,  0.25,  -0.25, // surf 3    3->4
        -0.25, 0.00,  -0.25, // surf 4    1->4
        0.00,  -0.25, 0.25,  // surf 5    5->6
        0.25,  0.00,  0.25,  // surf 6    6->7
        0.00,  0.25,  0.25,  // surf 7    7->8
        -0.25, 0.00,  0.25,  // surf 8    5->8
        -0.25, -0.25, 0.00,  // surf 9    1->5
        0.25,  -0.25, 0.00,  // surf 10   2->6
        0.25,  0.25,  0.00,  // surf 11   3->7
        -0.25, 0.25,  0.00}; // surf 12   4->8

    static constexpr double intgLocShift[numScsIp * dim] = {
        0.00,  -0.50, -0.50, // surf 1    1->2
        0.50,  0.00,  -0.50, // surf 2    2->3
        0.00,  0.50,  -0.50, // surf 3    3->4
        -0.50, 0.00,  -0.50, // surf 4    1->4
        0.00,  -0.50, 0.50,  // surf 5    5->6
        0.50,  0.00,  0.50,  // surf 6    6->7
        0.00,  0.50,  0.50,  // surf 7    7->8
        -0.50, 0.00,  0.50,  // surf 8    5->8
        -0.50, -0.50, 0.00,  // surf 9    1->5
        0.50,  -0.50, 0.00,  // surf 10   2->6
        0.50,  0.50,  0.00,  // surf 11   3->7
        -0.50, 0.50,  0.00}; // surf 12   4->8

    static KOKKOS_INLINE_FUNCTION void cofactorMatrix(scalar adjJac[][3],
                                                      const scalar jact[][3])
    {
        adjJac[0][0] = jact[1][1] * jact[2][2] - jact[2][1] * jact[1][2];
        adjJac[0][1] = jact[1][2] * jact[2][0] - jact[2][2] * jact[1][0];
        adjJac[0][2] = jact[1][0] * jact[2][1] - jact[2][0] * jact[1][1];

        adjJac[1][0] = jact[0][2] * jact[2][1] - jact[2][2] * jact[0][1];
        adjJac[1][1] = jact[0][0] * jact[2][2] - jact[2][0] * jact[0][2];
        adjJac[1][2] = jact[0][1] * jact[2][0] - jact[2][1] * jact[0][0];

        adjJac[2][0] = jact[0][1] * jact[1][2] - jact[1][1] * jact[0][2];
        adjJac[2][1] = jact[0][2] * jact[1][0] - jact[1][2] * jact[0][0];
        adjJac[2][2] = jact[0][0] * jact[1][1] - jact[1][0] * jact[0][1];
    }

    static KOKKOS_INLINE_FUNCTION void cofactorMatrix(scalar adjJac[][2],
                                                      const scalar jact[][2])
    {
        adjJac[0][0] = jact[1][1];
        adjJac[0][1] = -jact[1][0];
        adjJac[1][0] = -jact[0][1];
        adjJac[1][1] = jact[0][0];
    }

    template <typename ViewType>
    static KOKKOS_FUNCTION void subdivide_hex_8(ViewType coords,
                                                scalar coordv[27][3])
    {
        /**
         * Subdivide the coordinates of a hex8 element into 8 hexs along edge,
         * face, and volume midpoints
         */
        constexpr unsigned numBaseNodes = 8;

        for (unsigned n = 0; n < numBaseNodes; ++n)
        {
            coordv[n][0] = coords(n, 0);
            coordv[n][1] = coords(n, 1);
            coordv[n][2] = coords(n, 2);
        }

        // Face-by-face ordering for the subdivided hex.  This is different than
        // what is done for a multilinear Hex27 element, which has equivalent
        // nodal locations.
        for (unsigned d = 0; d < dim; ++d)
        {
            // Face 0
            coordv[8][d] = 0.5 * (coords(0, d) + coords(1, d));  // edge 1
            coordv[9][d] = 0.5 * (coords(1, d) + coords(2, d));  // edge 2
            coordv[10][d] = 0.5 * (coords(2, d) + coords(3, d)); // edge 3
            coordv[11][d] = 0.5 * (coords(3, d) + coords(0, d)); // edge 4

            coordv[12][d] = 0.25 * (coords(0, d) + coords(1, d) + coords(2, d) +
                                    coords(3, d));

            // Face 1
            coordv[13][d] = 0.5 * (coords(4, d) + coords(5, d)); // edge 5
            coordv[14][d] = 0.5 * (coords(5, d) + coords(6, d)); // edge 6
            coordv[15][d] = 0.5 * (coords(6, d) + coords(7, d)); // edge 7
            coordv[16][d] = 0.5 * (coords(7, d) + coords(4, d)); // edge 8

            coordv[17][d] = 0.25 * (coords(4, d) + coords(5, d) + coords(6, d) +
                                    coords(7, d));

            // Face 2
            coordv[18][d] = 0.5 * (coords(1, d) + coords(5, d)); // edge 9
            coordv[19][d] = 0.5 * (coords(0, d) + coords(4, d)); // edge 10

            coordv[20][d] = 0.25 * (coords(0, d) + coords(1, d) + coords(4, d) +
                                    coords(5, d));

            // Face 3
            coordv[21][d] = 0.5 * (coords(3, d) + coords(7, d)); // edge 11
            coordv[22][d] = 0.5 * (coords(2, d) + coords(6, d)); // edge 12

            coordv[23][d] = 0.25 * (coords(2, d) + coords(3, d) + coords(6, d) +
                                    coords(7, d));

            // Face 4
            coordv[24][d] = 0.25 * (coords(1, d) + coords(2, d) + coords(5, d) +
                                    coords(6, d));

            // Face 5
            coordv[25][d] = 0.25 * (coords(0, d) + coords(3, d) + coords(4, d) +
                                    coords(7, d));

            // Volume centroid
            coordv[26][d] = 0.;
            for (unsigned n = 0; n < numBaseNodes; ++n)
            {
                coordv[26][d] += coords(n, d);
            }
            coordv[26][d] *= 0.125;
        }
    }

    template <typename ViewType>
    static KOKKOS_FUNCTION void
    quad_area_by_triangulation(int ics,
                               const scalar areacoords[4][3],
                               ViewType& area)
    {
        /**
         * Form up the area vec consistently with the triangulation used
         * in the Grandy algorithm, on each subcontrol volume hex
         *
         * "Efficient computation of volume of
         * Hexahedral Cells", Jeffrey Grandy, LLNL, UCRL-ID-128886,
         *  October 30, 1997.
         */
        using ftype = typename ViewType::value_type;

        area(ics, 0) = 0.0;
        area(ics, 1) = 0.0;
        area(ics, 2) = 0.0;

        const ftype xmid[3] = {0.25 * (areacoords[0][0] + areacoords[1][0] +
                                       areacoords[2][0] + areacoords[3][0]),
                               0.25 * (areacoords[0][1] + areacoords[1][1] +
                                       areacoords[2][1] + areacoords[3][1]),
                               0.25 * (areacoords[0][2] + areacoords[1][2] +
                                       areacoords[2][2] + areacoords[3][2])};

        ftype r1[3] = {areacoords[0][0] - xmid[0],
                       areacoords[0][1] - xmid[1],
                       areacoords[0][2] - xmid[2]};
        for (int itriangle = 0; itriangle < 4; ++itriangle)
        {
            const int t_index = (itriangle + 1) % 4;
            const ftype r2[3] = {areacoords[t_index][0] - xmid[0],
                                 areacoords[t_index][1] - xmid[1],
                                 areacoords[t_index][2] - xmid[2]};

            area(ics, 0) += r1[1] * r2[2] - r2[1] * r1[2];
            area(ics, 1) += r1[2] * r2[0] - r2[2] * r1[0];
            area(ics, 2) += r1[0] * r2[1] - r2[0] * r1[1];

            r1[0] = r2[0];
            r1[1] = r2[1];
            r1[2] = r2[2];
        }
        area(ics, 0) *= 0.5;
        area(ics, 1) *= 0.5;
        area(ics, 2) *= 0.5;
    }

    // KOKKOS_INLINE_FUNCTION
    // scalar shapeFunctionSingle(const unsigned i, const unsigned j) const
    // {
    // }

    static void preCalcShapeFunction(SFViewHost& shape_fcn)
    {
        const scalar half = 0.50;
        const scalar one4th = 0.25;
        const scalar one8th = 0.125;
        constexpr unsigned npts = numScsIp;

        // TODO: ADD SWITCH FOR SHIFTED/UNSHIFTED!!!!
        // TODO: ADD SWITCH FOR SHIFTED/UNSHIFTED!!!!
        // TODO: ADD SWITCH FOR SHIFTED/UNSHIFTED!!!!
        // TODO: ADD SWITCH FOR SHIFTED/UNSHIFTED!!!!
        // TODO: ADD SWITCH FOR SHIFTED/UNSHIFTED!!!!
        // TODO: ADD SWITCH FOR SHIFTED/UNSHIFTED!!!!
        // TODO: ADD SWITCH FOR SHIFTED/UNSHIFTED!!!!

        for (unsigned j = 0; j < npts; ++j)
        {
            const scalar s1 = intgLocShift[j * dim + 0];
            const scalar s2 = intgLocShift[j * dim + 1];
            const scalar s3 = intgLocShift[j * dim + 2];

            shape_fcn(j, 0) = one8th + one4th * (-s1 - s2 - s3) +
                              half * (s2 * s3 + s3 * s1 + s1 * s2) -
                              s1 * s2 * s3;
            shape_fcn(j, 1) = one8th + one4th * (s1 - s2 - s3) +
                              half * (s2 * s3 - s3 * s1 - s1 * s2) +
                              s1 * s2 * s3;
            shape_fcn(j, 2) = one8th + one4th * (s1 + s2 - s3) +
                              half * (-s2 * s3 - s3 * s1 + s1 * s2) -
                              s1 * s2 * s3;
            shape_fcn(j, 3) = one8th + one4th * (-s1 + s2 - s3) +
                              half * (-s2 * s3 + s3 * s1 - s1 * s2) +
                              s1 * s2 * s3;
            shape_fcn(j, 4) = one8th + one4th * (-s1 - s2 + s3) +
                              half * (-s2 * s3 - s3 * s1 + s1 * s2) +
                              s1 * s2 * s3;
            shape_fcn(j, 5) = one8th + one4th * (s1 - s2 + s3) +
                              half * (-s2 * s3 + s3 * s1 - s1 * s2) -
                              s1 * s2 * s3;
            shape_fcn(j, 6) = one8th + one4th * (s1 + s2 + s3) +
                              half * (s2 * s3 + s3 * s1 + s1 * s2) +
                              s1 * s2 * s3;
            shape_fcn(j, 7) = one8th + one4th * (-s1 + s2 + s3) +
                              half * (s2 * s3 - s3 * s1 - s1 * s2) -
                              s1 * s2 * s3;
        }
    }

    static void preCalcDeriv(HexDerivViewHost& deriv)
    {
        constexpr unsigned npts = numScsIp;

        const double half = 1.0 / 2.0;
        const double one4th = 1.0 / 4.0;

        // TODO: ADD SWITCH FOR SHIFTED/UNSHIFTED!!!!
        // TODO: ADD SWITCH FOR SHIFTED/UNSHIFTED!!!!
        // TODO: ADD SWITCH FOR SHIFTED/UNSHIFTED!!!!
        // TODO: ADD SWITCH FOR SHIFTED/UNSHIFTED!!!!
        // TODO: ADD SWITCH FOR SHIFTED/UNSHIFTED!!!!
        // TODO: ADD SWITCH FOR SHIFTED/UNSHIFTED!!!!
        // TODO: ADD SWITCH FOR SHIFTED/UNSHIFTED!!!!

        for (unsigned j = 0; j < npts; ++j)
        {
            const double s1 = intgLocShift[j * dim + 0];
            const double s2 = intgLocShift[j * dim + 1];
            const double s3 = intgLocShift[j * dim + 2];

            const double s1s2 = s1 * s2;
            const double s2s3 = s2 * s3;
            const double s1s3 = s1 * s3;

            deriv(j, 0, 0) = half * (s3 + s2) - s2s3 - one4th;
            deriv(j, 1, 0) = half * (-s3 - s2) + s2s3 + one4th;
            deriv(j, 2, 0) = half * (-s3 + s2) - s2s3 + one4th;
            deriv(j, 3, 0) = half * (+s3 - s2) + s2s3 - one4th;
            deriv(j, 4, 0) = half * (-s3 + s2) + s2s3 - one4th;
            deriv(j, 5, 0) = half * (+s3 - s2) - s2s3 + one4th;
            deriv(j, 6, 0) = half * (+s3 + s2) + s2s3 + one4th;
            deriv(j, 7, 0) = half * (-s3 - s2) - s2s3 - one4th;

            deriv(j, 0, 1) = half * (s3 + s1) - s1s3 - one4th;
            deriv(j, 1, 1) = half * (s3 - s1) + s1s3 - one4th;
            deriv(j, 2, 1) = half * (-s3 + s1) - s1s3 + one4th;
            deriv(j, 3, 1) = half * (-s3 - s1) + s1s3 + one4th;
            deriv(j, 4, 1) = half * (-s3 + s1) + s1s3 - one4th;
            deriv(j, 5, 1) = half * (-s3 - s1) - s1s3 - one4th;
            deriv(j, 6, 1) = half * (s3 + s1) + s1s3 + one4th;
            deriv(j, 7, 1) = half * (s3 - s1) - s1s3 + one4th;

            deriv(j, 0, 2) = half * (s2 + s1) - s1s2 - one4th;
            deriv(j, 1, 2) = half * (s2 - s1) + s1s2 - one4th;
            deriv(j, 2, 2) = half * (-s2 - s1) - s1s2 - one4th;
            deriv(j, 3, 2) = half * (-s2 + s1) + s1s2 - one4th;
            deriv(j, 4, 2) = half * (-s2 - s1) + s1s2 + one4th;
            deriv(j, 5, 2) = half * (-s2 + s1) - s1s2 + one4th;
            deriv(j, 6, 2) = half * (s2 + s1) + s1s2 + one4th;
            deriv(j, 7, 2) = half * (s2 - s1) - s1s2 + one4th;
        }
    }

    template <typename CoordViewType>
    static KOKKOS_INLINE_FUNCTION void
    gradOpInline(const ConstHexDerivView& referenceGradWeights,
                 const CoordViewType& coords,
                 scalar weights[8][3],
                 const unsigned ip)
    {
        scalar jact[dim][dim];
        for (unsigned i = 0; i < dim; ++i)
            for (unsigned j = 0; j < dim; ++j)
                jact[i][j] = scalar(0.0);

        scalar refGrad[nodesPerElement][dim];
        for (unsigned n = 0; n < nodesPerElement; ++n)
        {
            for (unsigned i = 0; i < dim; ++i)
            {
                refGrad[n][i] = referenceGradWeights(ip, n, i);
            }
            for (unsigned i = 0; i < dim; ++i)
            {
                for (unsigned j = 0; j < dim; ++j)
                {
                    jact[i][j] += refGrad[n][j] * coords(n, i);
                }
            }
        }

        scalar adjJac[dim][dim];
        cofactorMatrix(adjJac, jact);

        scalar det = scalar(0.0);
        for (unsigned i = 0; i < dim; ++i)
            det += jact[i][0] * adjJac[i][0];
        // STK_NGP_ThrowAssertMsg(
        //     stk::simd::are_any(det > tiny_positive_value()),
        //     "Problem with Jacobian determinant");

        const scalar inv_detj = scalar(1.0) / det;

        for (unsigned n = 0; n < nodesPerElement; ++n)
        {
            for (unsigned i = 0; i < dim; ++i)
            {
                weights[n][i] = scalar(0.0);
                for (unsigned j = 0; j < dim; ++j)
                {
                    weights[n][i] += adjJac[i][j] * refGrad[n][j];
                }
                weights[n][i] *= inv_detj;
            }
        }
    }

    template<typename HexDerivViewType, typename NodeCoordsViewType>
    static KOKKOS_INLINE_FUNCTION
    void calcInvJac(const HexDerivViewType& derivs,
                           const NodeCoordsViewType& coords,
                           int ip,
                           scalar invJac[3][3])
    {
        // Compute Jacobian J = sum_n (dN_n/dξ_j * x_n^i)
        scalar jac[3][3] = {0.0};
        for (unsigned n = 0; n < 8; ++n) {
            for (unsigned i = 0; i < 3; ++i) {
                for (unsigned j = 0; j < 3; ++j) {
                    jac[i][j] += derivs(ip, n, j) * coords(n, i);
                }
            }
        }
    
        // Compute determinant
        scalar det =
            jac[0][0] * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1]) -
            jac[0][1] * (jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0]) +
            jac[0][2] * (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]);
    
        // Compute adjugate (cofactor matrix)
        scalar adj[3][3];
        adj[0][0] =  jac[1][1]*jac[2][2] - jac[1][2]*jac[2][1];
        adj[0][1] = -jac[0][1]*jac[2][2] + jac[0][2]*jac[2][1];
        adj[0][2] =  jac[0][1]*jac[1][2] - jac[0][2]*jac[1][1];
    
        adj[1][0] = -jac[1][0]*jac[2][2] + jac[1][2]*jac[2][0];
        adj[1][1] =  jac[0][0]*jac[2][2] - jac[0][2]*jac[2][0];
        adj[1][2] = -jac[0][0]*jac[1][2] + jac[0][2]*jac[1][0];
    
        adj[2][0] =  jac[1][0]*jac[2][1] - jac[1][1]*jac[2][0];
        adj[2][1] = -jac[0][0]*jac[2][1] + jac[0][1]*jac[2][0];
        adj[2][2] =  jac[0][0]*jac[1][1] - jac[0][1]*jac[1][0];
    
        // Inverse: invJac = adj / det
        scalar invDet = scalar(1.0) / det;
        for (unsigned i = 0; i < 3; ++i) {
            for (unsigned j = 0; j < 3; ++j) {
                invJac[i][j] = adj[i][j] * invDet;
            }
        }
    }

    template <typename CoordViewType>
    static KOKKOS_INLINE_FUNCTION auto
    calcInvJac(const ConstHexDerivView& referenceGradWeights,
               const CoordViewType& coords,
               const unsigned ip)
    {
#ifndef __CUDA_ARCH__
        Eigen::Matrix<scalar, nodesPerElement, dim> refGrad;
        for (unsigned n = 0; n < nodesPerElement; ++n)
        {
            for (unsigned d = 0; d < dim; ++d)
            {
                refGrad(n, d) = referenceGradWeights(ip, n, d);
            }
        }

        Eigen::Matrix<scalar, nodesPerElement, dim> nodeCoords;
        for (unsigned n = 0; n < nodesPerElement; ++n)
        {
            for (unsigned d = 0; d < dim; ++d)
            {
                nodeCoords(n, d) = coords(n, d);
            }
        }

        return (refGrad.transpose() * nodeCoords).inverse();

        // Eigen::Matrix3d jac = refGrad.transpose() * coords;
        // Eigen::Matrix3d invJac = jac.inverse();
        // return jac.inverse();
        // Eigen::Vector3d oneWeights =
        //     invJac * refGrad(1, Eigen::placeholders::all).transpose();
    #else
        return Eigen::Matrix3d::Zero();
    #endif
    }

#ifndef __CUDA_ARCH__

    template <typename CoordViewType>
    static KOKKOS_INLINE_FUNCTION void
    eigenGradOP(const ConstHexDerivView& referenceGradWeights,
                const CoordViewType& coords,
                Eigen::Matrix<scalar, 8, 3>& weights,
                const unsigned ip)
    {
        Eigen::Matrix<scalar, nodesPerElement, dim> refGrad;
        for (unsigned n = 0; n < nodesPerElement; ++n)
        {
            for (unsigned d = 0; d < dim; ++d)
            {
                refGrad(n, d) = referenceGradWeights(ip, n, d);
            }
        }

        Eigen::Matrix<scalar, nodesPerElement, dim> nodeCoords;
        for (unsigned n = 0; n < nodesPerElement; ++n)
        {
            for (unsigned d = 0; d < dim; ++d)
            {
                nodeCoords(n, d) = coords(n, d);
            }
        }

        Eigen::Matrix3d invJac = (refGrad.transpose() * nodeCoords).inverse();
        weights = (invJac * refGrad.transpose()).transpose();
    }
#endif

    template <typename CoordViewType, typename OutputViewType>
    static KOKKOS_INLINE_FUNCTION void
    gradOp(const ConstHexDerivView& referenceGradWeights,
           const CoordViewType& coords,
           OutputViewType& weights)
    {
        for (unsigned ip = 0; ip < referenceGradWeights.extent(0); ++ip)
        {
            scalar jact[dim][dim];
            for (unsigned i = 0; i < dim; ++i)
                for (unsigned j = 0; j < dim; ++j)
                    jact[i][j] = scalar(0.0);

            scalar refGrad[nodesPerElement][dim];
            for (unsigned n = 0; n < nodesPerElement; ++n)
            {
                for (unsigned i = 0; i < dim; ++i)
                {
                    refGrad[n][i] = referenceGradWeights(ip, n, i);
                }
                for (unsigned i = 0; i < dim; ++i)
                {
                    for (unsigned j = 0; j < dim; ++j)
                    {
                        jact[i][j] += refGrad[n][j] * coords(n, i);
                    }
                }
            }

            scalar adjJac[dim][dim];
            cofactorMatrix(adjJac, jact);

            scalar det = scalar(0.0);
            for (unsigned i = 0; i < dim; ++i)
                det += jact[i][0] * adjJac[i][0];
            // STK_NGP_ThrowAssertMsg(
            //     stk::simd::are_any(det > tiny_positive_value()),
            //     "Problem with Jacobian determinant");

            const scalar inv_detj = scalar(1.0) / det;

            for (unsigned n = 0; n < nodesPerElement; ++n)
            {
                for (unsigned i = 0; i < dim; ++i)
                {
                    weights(ip, n, i) = scalar(0.0);
                    for (unsigned j = 0; j < dim; ++j)
                    {
                        weights(ip, n, i) += adjJac[i][j] * refGrad[n][j];
                    }
                    weights(ip, n, i) *= inv_detj;
                }
            }
        }
    }

    template <typename ViewType>
    static KOKKOS_INLINE_FUNCTION void areaVectors(const ViewType& coords,
                                                   ViewType& areav)
    {
        constexpr int hex_edge_facet_table[12][4] = {{20, 8, 12, 26},
                                                     {24, 9, 12, 26},
                                                     {10, 12, 26, 23},
                                                     {11, 25, 26, 12},
                                                     {13, 20, 26, 17},
                                                     {17, 14, 24, 26},
                                                     {17, 15, 23, 26},
                                                     {16, 17, 26, 25},
                                                     {19, 20, 26, 25},
                                                     {20, 18, 24, 26},
                                                     {22, 23, 26, 24},
                                                     {21, 25, 26, 23}};

        scalar coordv[27][3];
        subdivide_hex_8(coords, coordv);

        constexpr unsigned npf = 4;
        constexpr unsigned nscs = numScsIp;
        for (unsigned ics = 0; ics < nscs; ++ics)
        {
            scalar scscoords[4][3];
            for (unsigned inode = 0; inode < npf; ++inode)
            {
                const int itrianglenode = hex_edge_facet_table[ics][inode];
                for (unsigned d = 0; d < dim; ++d)
                {
                    scscoords[inode][d] = coordv[itrianglenode][d];
                }
            }
            quad_area_by_triangulation(ics, scscoords, areav);
        }
    }
};

#endif