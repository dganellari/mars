#pragma once

#include <array>
#include <vector>
#include <cmath>

namespace mars
{
namespace fem
{

// Forward declarations
template<typename ElementTag, typename RealType>
class ReferenceElement;

// Quadrature point with weight
template<typename RealType>
struct QuadraturePoint
{
    RealType xi, eta, zeta; // Local coordinates
    RealType weight;
};

// Specialization for tetrahedra (linear Lagrange, P1)
template<typename RealType>
class ReferenceElement<TetTag, RealType>
{
public:
    static constexpr int numNodes      = 4;
    static constexpr int spatialDim    = 3;
    static constexpr int numQuadPoints = 4; // Degree 2 quadrature (exact for linear)

    // Evaluate basis function i at local coordinates (xi, eta, zeta)
    __host__ __device__ static RealType evaluateBasis(int i, RealType xi, RealType eta, RealType zeta)
    {
        switch (i)
        {
            case 0: return 1.0 - xi - eta - zeta; // N0
            case 1: return xi;                    // N1
            case 2: return eta;                   // N2
            case 3: return zeta;                  // N3
            default: return 0.0;
        }
    }

    // Evaluate gradient of basis function i (in reference space)
    // Returns (dN/dxi, dN/deta, dN/dzeta)
    __host__ __device__ static void evaluateGradient(int i, RealType grad[3])
    {
        switch (i)
        {
            case 0: // N0 = 1 - xi - eta - zeta
                grad[0] = -1.0;
                grad[1] = -1.0;
                grad[2] = -1.0;
                break;
            case 1: // N1 = xi
                grad[0] = 1.0;
                grad[1] = 0.0;
                grad[2] = 0.0;
                break;
            case 2: // N2 = eta
                grad[0] = 0.0;
                grad[1] = 1.0;
                grad[2] = 0.0;
                break;
            case 3: // N3 = zeta
                grad[0] = 0.0;
                grad[1] = 0.0;
                grad[2] = 1.0;
                break;
        }
    }

    // Get quadrature points and weights for tetrahedral element
    // Using degree 2 quadrature (4 points) - exact for linear elements
    static std::vector<QuadraturePoint<RealType>> getQuadraturePoints()
    {
        std::vector<QuadraturePoint<RealType>> qpts;

        // Symmetric 4-point rule for tetrahedra
        const RealType a = 0.585410196624968;
        const RealType b = 0.138196601125011;
        const RealType w = 1.0 / 24.0; // Weight (volume of ref tet = 1/6)

        qpts.push_back({b, b, b, w});
        qpts.push_back({a, b, b, w});
        qpts.push_back({b, a, b, w});
        qpts.push_back({b, b, a, w});

        return qpts;
    }

    // Device-compatible version: get individual quadrature point
    __host__ __device__ static QuadraturePoint<RealType> getQuadraturePoint(int q)
    {
        // Symmetric 4-point rule for tetrahedra
        const RealType a = 0.585410196624968;
        const RealType b = 0.138196601125011;
        const RealType w = 1.0 / 24.0; // Weight (volume of ref tet = 1/6)

        switch (q)
        {
            case 0: return {b, b, b, w};
            case 1: return {a, b, b, w};
            case 2: return {b, a, b, w};
            case 3: return {b, b, a, w};
            default: return {0, 0, 0, 0};
        }
    }

    // Compute Jacobian matrix for physical element
    // J[i][j] = dx_i/dxi_j where x = physical coords, xi = reference coords
    __host__ __device__ static void
    computeJacobian(const RealType* node_x, const RealType* node_y, const RealType* node_z, RealType J[3][3])
    {
        // For tetrahedra, Jacobian is constant (affine mapping)
        // J = [x1-x0, x2-x0, x3-x0]
        //     [y1-y0, y2-y0, y3-y0]
        //     [z1-z0, z2-z0, z3-z0]

        J[0][0] = node_x[1] - node_x[0];
        J[0][1] = node_x[2] - node_x[0];
        J[0][2] = node_x[3] - node_x[0];

        J[1][0] = node_y[1] - node_y[0];
        J[1][1] = node_y[2] - node_y[0];
        J[1][2] = node_y[3] - node_y[0];

        J[2][0] = node_z[1] - node_z[0];
        J[2][1] = node_z[2] - node_z[0];
        J[2][2] = node_z[3] - node_z[0];
    }

    // Compute determinant of Jacobian (6 * volume of tetrahedron)
    __host__ __device__ static RealType computeJacobianDeterminant(const RealType J[3][3])
    {
        return J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
               J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
    }

    // Compute inverse Jacobian (needed for gradient transformation)
    __host__ __device__ static void computeJacobianInverse(const RealType J[3][3], RealType Jinv[3][3])
    {
        RealType det    = computeJacobianDeterminant(J);
        RealType invDet = 1.0 / det;

        // Cofactor matrix / det
        Jinv[0][0] = invDet * (J[1][1] * J[2][2] - J[1][2] * J[2][1]);
        Jinv[0][1] = invDet * (J[0][2] * J[2][1] - J[0][1] * J[2][2]);
        Jinv[0][2] = invDet * (J[0][1] * J[1][2] - J[0][2] * J[1][1]);

        Jinv[1][0] = invDet * (J[1][2] * J[2][0] - J[1][0] * J[2][2]);
        Jinv[1][1] = invDet * (J[0][0] * J[2][2] - J[0][2] * J[2][0]);
        Jinv[1][2] = invDet * (J[0][2] * J[1][0] - J[0][0] * J[1][2]);

        Jinv[2][0] = invDet * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
        Jinv[2][1] = invDet * (J[0][1] * J[2][0] - J[0][0] * J[2][1]);
        Jinv[2][2] = invDet * (J[0][0] * J[1][1] - J[0][1] * J[1][0]);
    }

    // Transform gradients from reference to physical space
    // grad_phys = J^(-T) * grad_ref
    __host__ __device__ static void
    transformGradient(const RealType Jinv[3][3], const RealType grad_ref[3], RealType grad_phys[3])
    {
        grad_phys[0] = Jinv[0][0] * grad_ref[0] + Jinv[1][0] * grad_ref[1] + Jinv[2][0] * grad_ref[2];
        grad_phys[1] = Jinv[0][1] * grad_ref[0] + Jinv[1][1] * grad_ref[1] + Jinv[2][1] * grad_ref[2];
        grad_phys[2] = Jinv[0][2] * grad_ref[0] + Jinv[1][2] * grad_ref[1] + Jinv[2][2] * grad_ref[2];
    }
};

} // namespace fem
} // namespace mars
