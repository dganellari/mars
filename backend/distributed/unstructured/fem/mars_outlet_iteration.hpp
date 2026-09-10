#pragma once

#include <algorithm>
#include <cmath>

// dot = <R, delta R> and square = <delta R, delta R>, with fixed positive weights.
inline double outlet_correction_damping(double dot, double square, double maximum)
{
    if (!std::isfinite(dot) || !std::isfinite(square) || !std::isfinite(maximum)
        || !(dot < 0) || !(square > 0) || !(maximum > 0) || maximum > 1) return 0;
    return std::min(maximum, -dot / square);
}

inline bool outlet_correction_contracts(double before, double after, double damping)
{
    return std::isfinite(before) && std::isfinite(after) && damping > 0 && damping <= 1
        && after >= 0 && after < before && after <= before * (1.0 - 1e-4 * damping);
}
