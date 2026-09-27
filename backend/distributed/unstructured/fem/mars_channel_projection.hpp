#pragma once

#include <cmath>

inline bool channel_state_converged(double continuity, double balance, double change,
                                    double continuity_tolerance, double change_tolerance)
{
    return std::isfinite(continuity) && continuity >= 0 && continuity <= continuity_tolerance
        && std::isfinite(balance) && balance >= 0 && balance <= continuity_tolerance
        && std::isfinite(change) && change >= 0 && change <= change_tolerance;
}

template<class RealType>
__host__ __device__ inline bool channel_rectilinear_hex(const RealType coords[8][3], RealType tolerance)
{
    RealType low[3], high[3];
    for (int d = 0; d < 3; ++d) low[d] = high[d] = coords[0][d];
    for (int n = 0; n < 8; ++n)
        for (int d = 0; d < 3; ++d)
        {
            if (coords[n][d] < low[d]) low[d] = coords[n][d];
            if (coords[n][d] > high[d]) high[d] = coords[n][d];
        }
    for (int d = 0; d < 3; ++d) if (!(high[d] - low[d] > tolerance)) return false;
    unsigned corners = 0;
    for (int n = 0; n < 8; ++n)
    {
        int corner = 0;
        for (int d = 0; d < 3; ++d)
        {
            if (coords[n][d] - low[d] <= tolerance) continue;
            if (!(high[d] - coords[n][d] <= tolerance)) return false;
            corner |= 1 << d;
        }
        corners |= 1u << corner;
    }
    return corners == 255;
}

template<class RealType>
__host__ __device__ inline void channel_skew_flux(RealType flux, RealType q_left,
                                                 RealType q_right, RealType& left,
                                                 RealType& right)
{
    left = -RealType(0.5) * flux * q_right;
    right = RealType(0.5) * flux * q_left;
}

// The planar channel has no normal velocity DOFs. Fixed tangential DOFs are
// removed from both the pressure operator and its velocity correction.
template<class RealType>
__host__ __device__ inline void channel_project_gradient(bool fixed, RealType& gx,
                                                        RealType& gy, RealType& gz)
{
    if (fixed) { gx = RealType(0); gy = RealType(0); }
    gz = RealType(0);
}

template<class RealType>
__host__ __device__ inline RealType channel_lift_entry(bool fixed_row, bool fixed_column,
                                                      RealType value, RealType target)
{
    return !fixed_row && fixed_column ? -value * target : RealType(0);
}

template<class RealType>
__host__ __device__ inline RealType channel_opening_flux(RealType inlet_area,
                                                       RealType outlet_area,
                                                       RealType inlet_target,
                                                       RealType outlet_velocity)
{
    return inlet_area * inlet_target + outlet_area * outlet_velocity;
}
