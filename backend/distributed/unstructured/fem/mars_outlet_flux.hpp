#pragma once

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_OUTLET_HD __host__ __device__
#else
#define MARS_OUTLET_HD
#endif

// Kept independent of the CUDA runtime so gates execute the production arithmetic.
template<typename RealType>
MARS_OUTLET_HD inline void outlet_tet_gradient(const RealType coords[4][3], RealType& det,
                                              RealType dNdx[4][3])
{
    RealType J[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) J[i][j] = coords[j + 1][i] - coords[0][i];
    det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
        - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
        + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
    const RealType inv_det = RealType(1) / det;
    dNdx[1][0] =  (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * inv_det;
    dNdx[1][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]) * inv_det;
    dNdx[1][2] =  (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * inv_det;
    dNdx[2][0] = -(J[1][0] * J[2][2] - J[1][2] * J[2][0]) * inv_det;
    dNdx[2][1] =  (J[0][0] * J[2][2] - J[0][2] * J[2][0]) * inv_det;
    dNdx[2][2] = -(J[0][0] * J[1][2] - J[0][2] * J[1][0]) * inv_det;
    dNdx[3][0] =  (J[1][0] * J[2][1] - J[1][1] * J[2][0]) * inv_det;
    dNdx[3][1] = -(J[0][0] * J[2][1] - J[0][1] * J[2][0]) * inv_det;
    dNdx[3][2] =  (J[0][0] * J[1][1] - J[0][1] * J[1][0]) * inv_det;
    for (int i = 0; i < 3; ++i) dNdx[0][i] = -dNdx[1][i] - dNdx[2][i] - dNdx[3][i];
}

// The opposite node belongs to the gradient blend, not the face coefficient.
template<typename RealType>
MARS_OUTLET_HD inline RealType boundaryFaceCoefficient(RealType d0, RealType d1, RealType d2)
{
    return (d0 + d1 + d2) / RealType(3);
}

// One vertex sample carries A_f/3. The trace and reconstructed gradient are frozen.
template<typename RealType>
MARS_OUTLET_HD inline RealType outlet_sample_pressure_derivative(
    const RealType area[3], const RealType opposite_gradient[3], RealType coefficient)
{
    return -coefficient * (area[0] * opposite_gradient[0] + area[1] * opposite_gradient[1]
                           + area[2] * opposite_gradient[2]) / RealType(3);
}

template<typename RealType>
MARS_OUTLET_HD inline void outlet_facet_sample_flux(
    const RealType coords[4][3], const int face_nodes[3], int opposite,
    const RealType area[3], const RealType velocity[3][3], RealType p_opposite,
    const RealType trace[3], const RealType reconstructed_gradient[4][3],
    const RealType coefficient[3], RealType samples[3], RealType* sample_derivative = nullptr)
{
    RealType det, dNdx[4][3];
    outlet_tet_gradient(coords, det, dNdx);
    const RealType third = RealType(1) / RealType(3);
    const RealType d_face = boundaryFaceCoefficient(coefficient[0], coefficient[1], coefficient[2]);
    RealType stabilization = RealType(0);
    for (int i = 0; i < 3; ++i)
    {
        RealType compact_gradient = p_opposite * dNdx[opposite][i];
        RealType face_gradient = RealType(0);
        for (int r = 0; r < 3; ++r)
        {
            compact_gradient += trace[r] * dNdx[face_nodes[r]][i];
            face_gradient += reconstructed_gradient[face_nodes[r]][i];
        }
        const RealType blended_gradient =
            RealType(0.5) * (third * face_gradient + reconstructed_gradient[opposite][i]);
        stabilization += (blended_gradient - compact_gradient) * area[i];
    }
    stabilization *= d_face * third;
    for (int r = 0; r < 3; ++r)
        samples[r] = third * (velocity[r][0] * area[0] + velocity[r][1] * area[1]
                              + velocity[r][2] * area[2]) + stabilization;
    if (sample_derivative)
        *sample_derivative = outlet_sample_pressure_derivative(area, dNdx[opposite], d_face);
}

#undef MARS_OUTLET_HD
