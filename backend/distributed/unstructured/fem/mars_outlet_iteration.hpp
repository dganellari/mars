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

// For F = 0.5*rms(R)^2, slope = <R, delta R>_(1/V) / sum(V).
inline bool outlet_correction_contracts(double before, double after, double damping, double slope)
{
    if (!std::isfinite(before) || !std::isfinite(after) || !std::isfinite(damping)
        || !std::isfinite(slope) || !(before > 0) || !(slope < 0)
        || !(damping > 0) || damping > 1 || after < 0 || !(after < before)) return false;
    // Normalize without squaring a large/tiny norm; preserve a small measured decrease.
    const double relative_slope = (slope / before) / before;
    const double relative_change = ((after - before) / before) * (after / before + 1.0);
    return std::isfinite(relative_slope)
        && relative_change <= 2e-4 * damping * relative_slope;
}
