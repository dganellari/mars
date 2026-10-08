#pragma once
#include <cmath>
#include <limits>
#if defined(__CUDACC__)
#define MARS_COMPENSATED_HD __host__ __device__
#else
#define MARS_COMPENSATED_HD
#endif
namespace mars::segregated {
// Dot2 with an FMA product remainder preserves terms lost through cancellation.
// Explicit device rounding prevents contraction from removing the error terms.
struct CompensatedDot {
    double sum=0, error=0, low_magnitude=0, products=0;
    MARS_COMPENSATED_HD static double add_up(double a,double b) {
#if defined(__CUDA_ARCH__)
        return __dadd_ru(a,b);
#else
        return std::nextafter(a+b,std::numeric_limits<double>::infinity());
#endif
    }
    MARS_COMPENSATED_HD static double multiply_up(double a,double b) {
#if defined(__CUDA_ARCH__)
        return __dmul_ru(a,b);
#else
        return std::nextafter(a*b,std::numeric_limits<double>::infinity());
#endif
    }
    MARS_COMPENSATED_HD static double divide_up(double a,double b) {
#if defined(__CUDA_ARCH__)
        return __ddiv_ru(a,b);
#else
        return std::nextafter(a/b,std::numeric_limits<double>::infinity());
#endif
    }
    MARS_COMPENSATED_HD static double add(double a,double b) {
#if defined(__CUDA_ARCH__)
        return __dadd_rn(a,b);
#else
        return a+b;
#endif
    }
    MARS_COMPENSATED_HD void product(double a,double b) {
#if defined(__CUDA_ARCH__)
        const double p=__dmul_rn(a,b), remainder=__fma_rn(a,b,-p);
#else
        // The rounded product must exist separately even with FMA contraction.
        const volatile double rounded=a*b;
        const double p=rounded, remainder=std::fma(a,b,-p);
#endif
        const double next=add(sum,p), z=add(next,-sum);
        const double lost=add(add(sum,-add(next,-z)),add(p,-z));
        const double low=add(lost,remainder);
        low_magnitude=add(low_magnitude,fabs(low));
        products+=1;
        error=add(error,low);
        sum=next;
    }
    MARS_COMPENSATED_HD double value() const { return add(sum,error); }
    MARS_COMPENSATED_HD double error_bound() const {
        // Ogita-Rump-Oishi, Dot2Err (Algorithm 5.8), including underflow.
        // Count the exact leading zero product used by this streaming form.
        const double u=0.5*std::numeric_limits<double>::epsilon(), n=products+1;
        if (!(2*n*u<1)) return std::numeric_limits<double>::infinity();
        const double delta=divide_up(n*u,1-2*n*u);
        const double alpha=add_up(multiply_up(u,fabs(value())),
            add_up(multiply_up(delta,low_magnitude),3*(std::numeric_limits<double>::denorm_min()/u)));
        // Round upward instead of relying on the last division's rounding.
        return divide_up(alpha,1-2*u);
    }
};
} // namespace mars::segregated
#undef MARS_COMPENSATED_HD
