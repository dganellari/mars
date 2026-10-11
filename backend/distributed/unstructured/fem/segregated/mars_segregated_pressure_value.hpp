#pragma once
#include "mars_segregated_compensated_dot.hpp"
#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_PRESSURE_HD __host__ __device__
#else
#define MARS_PRESSURE_HD
#endif
namespace mars::segregated {
// Keep the update remainder until the pressure action has consumed it.
struct PressureValue {
    double high=0, low=0;
    MARS_PRESSURE_HD static PressureValue sum(double a,double b) {
        const double h=CompensatedDot::add(a,b), z=CompensatedDot::add(h,-a);
        return {h,CompensatedDot::add(CompensatedDot::add(a,-CompensatedDot::add(h,-z)),CompensatedDot::add(b,-z))};
    }
    MARS_PRESSURE_HD PressureValue operator+(PressureValue b) const {
        auto s=sum(high,b.high), t=sum(low,b.low);
        auto u=sum(s.low,t.high), v=sum(s.high,u.high);
        return sum(v.high,CompensatedDot::add(v.low,CompensatedDot::add(u.low,t.low)));
    }
    MARS_PRESSURE_HD PressureValue operator-() const { return {-high,-low}; }
    MARS_PRESSURE_HD PressureValue operator-(PressureValue b) const { return *this+(-b); }
    MARS_PRESSURE_HD PressureValue operator*(double b) const {
#if defined(__CUDA_ARCH__)
        const double h=__dmul_rn(high,b), l=__fma_rn(high,b,-h);
#else
        const volatile double product=high*b;
        const double h=product, l=std::fma(high,b,-h);
#endif
        return sum(h,CompensatedDot::add(l,low*b));
    }
    MARS_PRESSURE_HD PressureValue operator/(double b) const {
        const double q=high/b;
        const auto remainder=*this-PressureValue{q,0}*b;
        return sum(q,CompensatedDot::add(remainder.high,remainder.low)/b);
    }
    MARS_PRESSURE_HD double rounded() const { return CompensatedDot::add(high,low); }
    MARS_PRESSURE_HD static PressureValue load(const double* high,const double* low,int i) {
        return {high[i],low?low[i]:0};
    }
    MARS_PRESSURE_HD void store(double* h,double* l,int i) const { h[i]=high; l[i]=low; }
};
}
#undef MARS_PRESSURE_HD
