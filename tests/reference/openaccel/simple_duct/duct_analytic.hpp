#pragma once
// Fully developed laminar flow in a rectangular duct (Boussinesq series; White, Viscous Fluid
// Flow, eq. 3-48; Shah & London 1978). Cross-section |y|<=W/2, |z|<=H/2, flow along +x,
// G=-dp/dx>0:
//
//   mu (u_yy + u_zz) = -G,  u = 0 on the walls.
//
// With a the shorter and b the longer half-side, s the coordinate across the short side and
// t the one along the long side, lambda_i = i pi/(2a), i odd:
//
//   u = G/(2mu) (a^2 - s^2)
//       - 16 a^2 G/(mu pi^3) sum_i (-1)^((i-1)/2) i^-3 cos(lambda_i s) cosh(lambda_i t)/cosh(lambda_i b)
//
//   U_mean = a^2 G/(3mu) K,   K = 1 - 192 a/(pi^5 b) sum_i tanh(i pi b/(2a)) / i^5
//
// The first term is the plane-Poiseuille profile across the short side; the series corrects it
// near the short walls |t|=b, and its terms decay like exp(-lambda_i (b-|t|)). The same solution
// with the roles of the sides exchanged (a<->b, s<->t) decays like exp(-i pi (a-|s|)/(2b)), so
// each point uses the faster of the two. cosh(lambda t)/cosh(lambda b) is evaluated as
// exp(lambda(|t|-b)) (1+exp(-2 lambda |t|))/(1+exp(-2 lambda b)), which never overflows. The
// sum stops when the remaining terms are below 1e-18 of the leading scale, or after 40000
// terms (tail below 1/(4 N^2) ~ 2e-10 of it; reached only within ~1e-4 a of a corner). Points on
// a wall return 0. Host-only, no MARS dependencies.
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace duct {

constexpr double pi=3.14159265358979323846;

struct Duct {
    double width=2, height=1;   // along y and z
    double viscosity=.1;        // mu

    double a() const { return .5*std::min(width,height); }   // shorter half-side
    double b() const { return .5*std::max(width,height); }   // longer half-side
    double area() const { return width*height; }
    double hydraulic_diameter() const { return 2*width*height/(width+height); }   // 4A/P

    // K(b/a) of U_mean = a^2 G K/(3 mu); 1 for parallel plates, 0.4217 for the square.
    double shape_factor() const {
        const double ratio=b()/a(); double sum=0;
        for (int i=1;i<200000;i+=2) {
            const double term=std::tanh(i*pi*ratio/2)/std::pow(double(i),5);
            sum+=term; if (term<1e-19*sum) break;
        }
        return 1-192/(std::pow(pi,5)*ratio)*sum;
    }
    // G that carries the mean velocity U, and its inverse.
    double pressure_gradient(double mean_velocity) const { return 3*viscosity*mean_velocity/(a()*a()*shape_factor()); }
    double mean_velocity(double G) const { return a()*a()*G*shape_factor()/(3*viscosity); }
    double flow_rate(double G) const { return mean_velocity(G)*area(); }
    // Fanning f Re_Dh = G D_h^2 / (2 mu U): 14.227 square, 15.548 for aspect 1:2, 24 for plates.
    double fanning_re() const { const double d=hydraulic_diameter(); return 3*d*d/(2*a()*a()*shape_factor()); }

    // u(y,z) for pressure gradient G.
    double velocity(double y,double z,double G) const {
        const bool z_short=height<=width;
        const double s=std::abs(z_short?z:y), t=std::abs(z_short?y:z), A=a(), B=b();
        if (s>=A || t>=B) return 0;
        return (B-t)/A>=(A-s)/B ? expansion(s,t,A,B,G) : expansion(t,s,B,A,G);
    }
    // Plane profile across the half-width `half_across`, plus the series decaying along `along`.
    double expansion(double across,double along,double half_across,double half_along,double G) const {
        double series=0;
        for (int i=1;i<80000;i+=2) {
            const double lambda=i*pi/(2*half_across);
            const double decay=std::exp(lambda*(along-half_along));
            const double ratio=decay*(1+std::exp(-2*lambda*along))/(1+std::exp(-2*lambda*half_along));
            const double sign=((i-1)/2)%2==0?1:-1;
            series+=sign*std::cos(lambda*across)*ratio/(double(i)*i*i);
            if (2*decay/(double(i)*i*i)<1e-18) break;
        }
        return G/viscosity*(.5*(half_across*half_across-across*across)-16*half_across*half_across/(pi*pi*pi)*series);
    }
    double centerline(double G) const { return velocity(0,0,G); }
};

// Durst et al. (J. Fluids Eng. 127, 2005) pipe entrance length, L/D = [0.619^1.6 + (0.0567 Re)^1.6]^(1/1.6),
// applied with D = D_h: an order-of-magnitude estimate for a rectangular duct. The comparator
// (duct_compare.py) measures the development length of the discrete solution instead.
inline double durst_entrance_length(double reynolds,double diameter) {
    return diameter*std::pow(std::pow(.619,1.6)+std::pow(.0567*reynolds,1.6),1/1.6);
}

} // namespace duct
