#pragma once
#include "mars_segregated_simple.hpp"
#include <algorithm>
#include <cmath>

#if defined(__CUDACC__)
#define MARS_METRICS_HD __host__ __device__
#else
#define MARS_METRICS_HD
#endif
namespace mars::segregated {
struct SimpleSums {
    double volume=0, momentum2=0, continuity2=0, continuity=0;
    double velocity_change2=0, pressure_change2=0, speed2=0;
    double inlet=0, outlet=0, inlet_area=0, flux_change=0;
    int closed=0, changed=0, invalid=0;
};
struct SimpleSumCombine {
    MARS_METRICS_HD SimpleSums operator()(SimpleSums a, const SimpleSums& b) const {
        a.volume+=b.volume; a.momentum2+=b.momentum2; a.continuity2+=b.continuity2; a.continuity+=b.continuity;
        a.velocity_change2+=b.velocity_change2; a.pressure_change2+=b.pressure_change2;
        if (b.speed2>a.speed2) a.speed2=b.speed2;
        if (b.flux_change>a.flux_change) a.flux_change=b.flux_change;
        a.inlet+=b.inlet; a.outlet+=b.outlet; a.inlet_area+=b.inlet_area;
        a.closed+=b.closed; a.changed+=b.changed; a.invalid+=b.invalid; return a;
    }
};
struct SimpleNodeSums {
    SimpleState state; const double *rhs,*old_velocity,*old_pressure;
    MARS_METRICS_HD SimpleSums operator()(int n) const {
        SimpleSums a; const double v=state.volume[n];
        if (!(v>0) || !finite_coefficient(v)) { a.invalid=1; return a; }
        a.volume=v; a.continuity=state.mass_divergence[n];
        a.continuity2=a.continuity*a.continuity/v;
        // Undo boundary RHS relaxation to measure the physical momentum defect.
        const double factor=state.boundary_factor[n]>0?.75:1.;
        for (int j=0;j<3;++j) {
            const double r=rhs[3*n+j]/factor, u=state.velocity[3*n+j];
            const double du=u-old_velocity[3*n+j];
            a.momentum2+=r*r/v; a.velocity_change2+=v*du*du; a.speed2+=u*u;
        }
        const double dp=state.pressure[n]-old_pressure[n]; a.pressure_change2=v*dp*dp;
        a.invalid=!finite_coefficient(a.momentum2) || !finite_coefficient(a.continuity2)
            || !finite_coefficient(a.velocity_change2) || !finite_coefficient(a.pressure_change2)
            || !finite_coefficient(a.speed2);
        return a;
    }
};
struct SimpleFaceSums {
    SimpleMesh mesh; SimpleState state; const int* old_flags;
    MARS_METRICS_HD SimpleSums operator()(int i) const {
        SimpleSums a; const auto f=mesh.faces[i];
        for (int j=0;j<3;++j) {
            const double q=state.boundary_flux[3*i+j];
            a.invalid+=!finite_coefficient(q);
            if (f.kind==0) a.inlet+=q;
            else if (f.kind==1) a.outlet+=q;
            else a.invalid+=q!=0;
            a.invalid+=state.reversal[3*i+j]!=state.reversal[3*i];
        }
        if (f.kind==0) {
            double area[3]; tet_boundary_area(mesh.geometry[f.element],f.ordinal,area);
            a.inlet_area=3*sqrt(area[0]*area[0]+area[1]*area[1]+area[2]*area[2]);
        }
        if (f.kind==1) { a.closed=state.reversal[3*i]!=0; a.changed=state.reversal[3*i]!=old_flags[3*i]; }
        return a;
    }
};
struct SimpleFluxChange {
    const double *flux,*old_flux;
    MARS_METRICS_HD SimpleSums operator()(int i) const {
        SimpleSums a; const double dq=flux[i]-old_flux[i];
        a.flux_change=dq<0?-dq:dq; a.invalid=!finite_coefficient(dq); return a;
    }
};
struct SimpleMetrics {
    double momentum, continuity, flux, velocity_change, pressure_change, cancellation;
    bool finite;
    double flux_change=0;
};
// L=1 m is the public channel height. These are MARS norms, not OpenAccel's RMS columns.
inline SimpleMetrics simple_metrics(const SimpleSums& s,const SimpleControls& c,double length=1) {
    const double target=c.density*c.inlet_speed*s.inlet_area;
    SimpleMetrics m{sqrt(s.momentum2/s.volume)*length/(c.density*c.inlet_speed*c.inlet_speed),
        sqrt(s.continuity2/s.volume)*length/(c.density*c.inlet_speed),
        std::abs(s.inlet+s.outlet)/target, sqrt(s.velocity_change2/s.volume)/c.inlet_speed,
        sqrt(s.pressure_change2/s.volume)/(c.density*c.inlet_speed*c.inlet_speed),
        std::abs(s.continuity-s.inlet-s.outlet)/target, false, s.flux_change/target};
    m.finite=!s.invalid && s.volume>0 && target>0 && length>0;
    for (double value:{m.momentum,m.continuity,m.flux,m.velocity_change,m.pressure_change,m.cancellation,m.flux_change})
        m.finite=m.finite && std::isfinite(value);
    return m;
}
inline bool simple_converged(const SimpleMetrics& m,int iterations,int changed,
    double residual_tolerance,double flux_tolerance,double change_tolerance) {
    return m.finite && iterations>=2 && changed==0 && m.momentum<=residual_tolerance
        && m.continuity<=residual_tolerance && m.flux<=flux_tolerance
        && m.velocity_change<=change_tolerance && m.pressure_change<=change_tolerance
        && m.flux_change<=change_tolerance && m.cancellation<=1e-10;
}
}
#undef MARS_METRICS_HD
