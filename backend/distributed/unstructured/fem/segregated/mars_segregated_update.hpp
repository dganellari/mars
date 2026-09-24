#pragma once
#include <cmath>
#include <limits>

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_UPDATE_HD __host__ __device__
#else
#define MARS_UPDATE_HD
#endif

namespace mars::segregated {

MARS_UPDATE_HD inline double pressure_update(double old, double increment, double alpha)
{
    return old + alpha*increment;
}

// d already contains momentum relaxation; grad_increment is not pressure-relaxed.
MARS_UPDATE_HD inline void velocity_update(const double* old, const double* d,
                                           const double* grad_increment, double* result)
{
    for (int j = 0; j < 3; ++j) result[j] = old[j] - d[j]*grad_increment[j];
}

struct FluxUpdate {
    double density, velocity[3], influence[3], compact_gradient[3], reconstructed_gradient[3];
    double original_force[3], reconstructed_force[3], area[3], old, alpha;
};

// Fixed-frame mass flux in kg/s. Interpolation and gradient reconstruction are inputs.
MARS_UPDATE_HD inline double mass_flux_update(const FluxUpdate& x, bool outlet, bool reversed)
{
    if (outlet && reversed) return 0;
    double fresh = 0;
    for (int j = 0; j < 3; ++j) {
        if (outlet) {
            fresh += (x.density*x.velocity[j] - x.density*x.influence[j]
                      *(x.compact_gradient[j]-x.reconstructed_gradient[j]))*x.area[j];
        } else {
            fresh += x.density*x.velocity[j]*x.area[j];
            fresh -= x.density*x.influence[j]*(x.compact_gradient[j]-x.reconstructed_gradient[j])*x.area[j];
        }
        fresh += x.density*x.influence[j]*(x.original_force[j]-x.reconstructed_force[j])*x.area[j];
    }
    return x.alpha*fresh + (1-x.alpha)*x.old;
}

MARS_UPDATE_HD inline double inlet_flux_update(double density, const double* prescribed,
                                               const double* area, double old, double alpha)
{
    double fresh = 0;
    for (int j = 0; j < 3; ++j) fresh += (density*prescribed[j])*area[j];
    return alpha*fresh + (1-alpha)*old;
}

// Flags describe one artificial wall per Tri3 face, not three independent valves.
// The caller supplies a nonzero face normal and uniform face flags.
MARS_UPDATE_HD inline void outlet_reversal_update(const double* old_flux, const int* old_flags,
                                                  const double* velocity, const double* pressure,
                                                  const double* trace, const double* area, bool ignore,
                                                  double* flux, int* flags)
{
    double net = 0, normal[3]{}, mean_velocity[3]{}, area_squared = 0;
    for (int s = 0; s < 3; ++s) {
        net += old_flux[s];
        flux[s] = old_flux[s] > 0 ? old_flux[s] : 0;
        flags[s] = old_flags[s];
    }
    for (int j = 0; j < 3; ++j) {
        for (int s = 0; s < 3; ++s) normal[j] += area[3*s+j];
        area_squared += normal[j]*normal[j];
    }
    const double magnitude = ::sqrt(area_squared);
    double mean_pressure = 0, mean_trace = 0, normal_velocity = 0;
    for (int s = 0; s < 3; ++s) {
        mean_pressure += (1.0/3)*pressure[s];
        for (int j = 0; j < 3; ++j) mean_velocity[j] += (1.0/3)*velocity[3*s+j];
        mean_trace += trace[s];
    }
    mean_trace /= 3;
    for (int j = 0; j < 3; ++j) normal_velocity += mean_velocity[j]*normal[j]/magnitude;
    if (old_flags[0] == 0) {
        if (net < -std::numeric_limits<double>::epsilon())
            for (int s = 0; s < 3; ++s) { if (!ignore) flags[s] = 1; flux[s] = 0; }
    } else if (normal_velocity >= 0 && mean_trace <= mean_pressure) {
        if (!ignore) for (int s = 0; s < 3; ++s) flags[s] = 0;
    } else {
        for (int s = 0; s < 3; ++s) flux[s] = 0;
    }
}

MARS_UPDATE_HD inline void outlet_trace_moment(const double* pressure, const double* shape,
                                              const double* area, bool reversed, double* moment)
{
    moment[0] = moment[1] = 0;
    if (reversed) return;
    double squared = 0, sample = 0;
    for (int j = 0; j < 3; ++j) squared += area[j]*area[j];
    for (int s = 0; s < 3; ++s) sample += shape[s]*pressure[s];
    moment[1] = ::sqrt(squared);
    moment[0] = sample*moment[1];
}

MARS_UPDATE_HD inline double outlet_trace_update(double nearest, double prescribed, double mean,
                                                 double beta, double old, bool reversed)
{
    return reversed ? old : prescribed + (1-beta)*(nearest-mean);
}
} // namespace mars::segregated
#undef MARS_UPDATE_HD
