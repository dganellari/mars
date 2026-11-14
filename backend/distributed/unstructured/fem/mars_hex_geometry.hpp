#pragma once

// Geometry of the trilinear (Q1) hex: the Jacobian determinant and the volumes
// of the eight sub-control volumes (the lumped mass of CVFEM). Shared by the
// hex Navier-Stokes solver (mars_navier_stokes.hpp) and the pump solvers
// (mars_ns_pump_solver.hpp, mars_ns_solver.hpp).
//
// Corner order (Exodus HEX8): (-,-,-) (+,-,-) (+,+,-) (-,+,-) (-,-,+) (+,-,+) (+,+,+) (-,+,+).

#include "../utils/mars_gpu_runtime.hpp"

#include <cmath>

namespace mars
{
namespace fem
{

// det J of the trilinear map of a hex at reference point xi in [-1, 1]^3.
template<typename RealType>
__host__ __device__ inline RealType hexJacobianDeterminant(const RealType x[8][3], const RealType xi[3])
{
    const int sign[8][3] = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
                            {-1, -1, 1},  {1, -1, 1},  {1, 1, 1},  {-1, 1, 1}};
    RealType J[3][3] = {};
    for (int n = 0; n < 8; ++n)
    {
        RealType f[3] = {1 + sign[n][0] * xi[0], 1 + sign[n][1] * xi[1], 1 + sign[n][2] * xi[2]};
        RealType dN[3] = {sign[n][0] * f[1] * f[2] / 8, sign[n][1] * f[0] * f[2] / 8, sign[n][2] * f[0] * f[1] / 8};
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                J[a][b] += x[n][a] * dN[b];
    }
    return J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
           J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
}

// Volume of the sub-control volume of corner c: the image of the reference
// octant between the corner and the element centre. det J is quadratic in each
// reference coordinate, so the 2x2x2 Gauss rule on the octant is exact.
template<typename RealType>
__host__ __device__ inline RealType hexSubVolume(const RealType x[8][3], int c)
{
    const int sign[8][3] = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
                            {-1, -1, 1},  {1, -1, 1},  {1, 1, 1},  {-1, 1, 1}};
    const RealType g[2]  = {RealType(0.5) - RealType(0.5) / sqrt(RealType(3)),
                            RealType(0.5) + RealType(0.5) / sqrt(RealType(3))};
    RealType volume      = 0;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k)
            {
                RealType xi[3] = {sign[c][0] * g[i], sign[c][1] * g[j], sign[c][2] * g[k]};
                volume += hexJacobianDeterminant(x, xi);
            }
    return volume / 8; // octant volume 1, eight points
}

// The trilinear basis N and its physical gradients dNdx at reference point xi
// in [-1, 1]^3; returns det J (<= 0 for an inverted or degenerate point, where
// dNdx is left unset).
template<typename RealType>
__host__ __device__ inline RealType hexBasisGradients(const RealType x[8][3], const RealType xi[3], RealType N[8],
                                                      RealType dNdx[8][3])
{
    const int sign[8][3] = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
                            {-1, -1, 1},  {1, -1, 1},  {1, 1, 1},  {-1, 1, 1}};
    RealType dN[8][3];
    RealType J[3][3] = {};
    for (int n = 0; n < 8; ++n)
    {
        RealType f[3] = {1 + sign[n][0] * xi[0], 1 + sign[n][1] * xi[1], 1 + sign[n][2] * xi[2]};
        N[n]     = f[0] * f[1] * f[2] / 8;
        dN[n][0] = sign[n][0] * f[1] * f[2] / 8;
        dN[n][1] = sign[n][1] * f[0] * f[2] / 8;
        dN[n][2] = sign[n][2] * f[0] * f[1] / 8;
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                J[a][b] += x[n][a] * dN[n][b];
    }
    const RealType det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) -
                         J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
                         J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
    if (!(det > RealType(0))) return det;
    // Jinv[b][a] = d xi_b / d x_a
    RealType Jinv[3][3];
    Jinv[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) / det;
    Jinv[0][1] = (J[0][2] * J[2][1] - J[0][1] * J[2][2]) / det;
    Jinv[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) / det;
    Jinv[1][0] = (J[1][2] * J[2][0] - J[1][0] * J[2][2]) / det;
    Jinv[1][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) / det;
    Jinv[1][2] = (J[0][2] * J[1][0] - J[0][0] * J[1][2]) / det;
    Jinv[2][0] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) / det;
    Jinv[2][1] = (J[0][1] * J[2][0] - J[0][0] * J[2][1]) / det;
    Jinv[2][2] = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) / det;
    for (int n = 0; n < 8; ++n)
        for (int a = 0; a < 3; ++a)
            dNdx[n][a] = dN[n][0] * Jinv[0][a] + dN[n][1] * Jinv[1][a] + dN[n][2] * Jinv[2][a];
    return det;
}

} // namespace fem
} // namespace mars
