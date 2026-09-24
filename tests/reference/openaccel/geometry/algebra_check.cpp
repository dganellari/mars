#include "mars_segregated_geometry.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace mars::segregated;
int checks = 0;
void close(double actual, double expected) {
    ++checks;
    if (!std::isfinite(actual) || std::abs(actual-expected) > 2e-12*std::max(1.,std::abs(expected)))
        throw std::runtime_error("geometry algebra mismatch");
}
void require(bool condition) { ++checks; if (!condition) throw std::runtime_error("geometry validity mismatch"); }

void affine_fixture(const double* x)
{
    TetGeometry<double> g; require(tet_geometry(x,g));
    double field[4], constant[4]{2,2,2,2}, numerator[12], closure[12]{};
    const double exact[]{2,-3,5};
    for (int n = 0; n < 4; ++n) field[n] = 1+2*x[3*n]-3*x[3*n+1]+5*x[3*n+2];
    for (int j = 0; j < 3; ++j) {
        double derivative = 0, sum = 0;
        for (int n = 0; n < 4; ++n) { derivative += field[n]*g.gradient[3*n+j]; sum += g.gradient[3*n+j]; }
        close(derivative,exact[j]); close(sum,0);
    }
    for (int s = 0; s < 6; ++s) {
        int l=tet_edge_node(s,0), r=tet_edge_node(s,1); double dot=0;
        for (int j = 0; j < 3; ++j) {
            closure[3*l+j] += g.area[3*s+j]; closure[3*r+j] -= g.area[3*s+j];
            dot += g.area[3*s+j]*(x[3*r+j]-x[3*l+j]);
        }
        close(dot,g.volume/2);
        double shape[4]; tet_sample_shape(s,false,shape);
        for (int n=0;n<4;++n) close(shape[n],(n==l||n==r)?13./36:5./36);
        tet_sample_shape(s,true,shape);
        for (int n=0;n<4;++n) close(shape[n],(n==l||n==r)?.5:0);
    }
    tet_gradient_numerator<1>(g,field,false,true,numerator);
    for (int f=0;f<4;++f) {
        double area[3], values[3], local[9]; tet_boundary_area(g,f,area);
        const int a=tet_face_node(f,0), b=tet_face_node(f,1), c=tet_face_node(f,2);
        for (int j=0;j<3;++j) {
            int k=(j+1)%3, l=(j+2)%3;
            close(area[j],((x[3*b+k]-x[3*a+k])*(x[3*c+l]-x[3*a+l])
                         -(x[3*b+l]-x[3*a+l])*(x[3*c+k]-x[3*a+k]))/6);
        }
        for (int n=0;n<3;++n) values[n]=field[tet_face_node(f,n)];
        tri_gradient_numerator<1>(area,values,false,true,local);
        for (int n=0;n<3;++n) for (int j=0;j<3;++j) {
            int index=3*tet_face_node(f,n)+j;
            closure[index]+=area[j]; numerator[index]+=local[3*n+j];
        }
        tri_gradient_numerator<1>(area,values,true,true,local);
        for (double value:local) close(value,0);
    }
    for (int n=0;n<4;++n) for (int j=0;j<3;++j) {
        close(closure[3*n+j],0); close(numerator[3*n+j]/(g.volume/4),exact[j]);
    }
    tet_gradient_numerator<1>(g,constant,false,true,numerator);
    for (double value:numerator) close(value,0);
    // Component layout must reproduce three independent scalar reconstructions.
    double vector_field[12], vector_num[36];
    for (int n=0;n<4;++n) for (int c=0;c<3;++c) vector_field[3*n+c]=(c+1)*field[n];
    tet_gradient_numerator<3>(g,vector_field,true,true,vector_num);
    tet_gradient_numerator<1>(g,field,true,true,numerator);
    for (int n=0;n<4;++n) for (int c=0;c<3;++c) for (int j=0;j<3;++j)
        close(vector_num[9*n+3*c+j],(c+1)*numerator[3*n+j]);
}

int main() {
    try {
        const double unit[]{0,0,0,1,0,0,0,1,0,0,0,1};
        affine_fixture(unit);
        TetGeometry<double> g; require(tet_geometry(unit,g)); close(g.volume,1./6);
        const double area[]{1./12,1./24,1./24};
        for (int j=0;j<3;++j) close(g.area[j],area[j]);
        double field[]{0,2,-3,5}, numerator[12];
        tet_gradient_numerator<1>(g,field,true,true,numerator);
        close(numerator[0]/(g.volume/4),3); close(numerator[1]/(g.volume/4),.5);
        close(numerator[2]/(g.volume/4),4.5); // shifted boundary reconstruction is not affine-exact
        for (int sample=0;sample<24;++sample) {
            double x[12]; const double scale=.25+.1*sample;
            for (int n=0;n<4;++n) {
                x[3*n] = 2+scale*(unit[3*n]+.2*unit[3*n+1]-.3*unit[3*n+2]);
                x[3*n+1] = -3+scale*(2*unit[3*n+1]+.4*unit[3*n+2]);
                x[3*n+2] = 1+scale*.7*unit[3*n+2];
            }
            affine_fixture(x);
        }
        double invalid[12]; std::copy(unit,unit+12,invalid);
        std::swap_ranges(invalid+3,invalid+6,invalid+6); require(!tet_geometry(invalid,g));
        std::copy(unit,unit+12,invalid); invalid[11]=0; require(!tet_geometry(invalid,g));
        invalid[0]=std::numeric_limits<double>::infinity(); require(!tet_geometry(invalid,g));
        double result=0, nan=std::numeric_limits<double>::quiet_NaN();
        require(finish_gradient(6.,2.,nan,.3,false,result)); close(result,3);
        require(finish_gradient(6.,2.,nan,1.,true,result)); close(result,3);
        require(finish_gradient(6.,2.,1.,.25,true,result)); close(result,1.5);
        require(!finish_gradient(1.,0.,0.,1.,false,result));
        require(!finish_gradient(nan,1.,0.,1.,false,result));
        require(!finish_gradient(1.,1.,0.,2.,false,result));
        std::cout << "PASS: " << checks << " independent geometry/reconstruction algebra checks\n";
        return 0;
    } catch (const std::exception& error) { std::cerr << "FAIL: " << error.what() << '\n'; return 1; }
}
