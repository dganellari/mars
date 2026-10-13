#pragma once
#include "mars_segregated_simple.hpp"

// Exercise the production boundary input on either host or device.
struct WallDistanceCheck {
    int* passed;
#if defined(__CUDACC__)
    __host__ __device__
#endif
    void operator()(int) const {
        using namespace mars::segregated;
        const double small=std::numeric_limits<double>::epsilon();
        const double scales[]={1.,2*small,small,.5*small};
        bool ok=true;
        for (double h:scales) {
            const double xyz[]={0,0,0, h,0,0, 0,h,0, 0,0,h};
            TetGeometry<double> g;
            ok=tet_geometry(xyz,g) && ok;
            const int nodes[]={0,1,2,3};
            double u[12]{},p[4]{},trace[3]{},flux[3]{};
            SimpleControls c;
            for (int face=0;face<4;++face) {
                double area[3]; tet_boundary_area(g,face,area);
                double norm=0,distance=0;
                const int opposite=tet_opposite_node(face),first=tet_face_node(face,0);
                for (int j=0;j<3;++j) {
                    norm+=area[j]*area[j];
                    distance+=(xyz[3*opposite+j]-xyz[3*first+j])*area[j];
                }
                norm=sqrt(norm); distance=fabs(distance)/norm;
                const double expected=c.viscosity*norm/(.25*(distance<small?small:distance));
                for (int initialized=0;initialized<2;++initialized) {
                    const auto boundary=simple_boundary(true,{0,face,2},nodes,g,u,p,trace,flux,c,initialized,nullptr);
                    for (double coefficient:boundary.values.wall_coefficient)
                        ok=ok && (initialized?fabs(coefficient/expected-1)<2e-14:coefficient==0);
                }
            }
        }
        *passed=ok;
    }
};
