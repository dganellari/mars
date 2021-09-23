#ifndef MANDELBROT_HPP
#define MANDELBROT_HPP

#include <complex>
#include "Image.hpp"

namespace mars {
    class Mandelbrot {
    public:
        Mandelbrot(){};
        int maxiter{200};
        int outofbounds{3};

        int operator()(double cx, double cy, double cz) const {
            std::complex<double> c{cx, cy};
            std::complex<double> z = c;
            int i = 0;
            while (i < maxiter) {
                double n = std::norm(z);
                if (n > outofbounds) {
                    break;
                }

                z = z * z + c;
                i++;
            }

            return i;
        }
    };
}  // namespace mars

#endif