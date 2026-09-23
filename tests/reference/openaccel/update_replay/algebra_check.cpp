#include "mars_segregated_update.hpp"
#include <algorithm>
#include <iostream>
#include <stdexcept>
using namespace mars::segregated;
int checks=0;
void eq(double a,double b) {
    ++checks;
    if(!std::isfinite(a) || std::abs(a-b)>2e-13*std::max(1.0,std::abs(b)))
        throw std::runtime_error("independent update algebra mismatch");
}
int main() {
    try {
        eq(pressure_update(10,8,.3),12.4);
        double u[]={1,2,3},d[]={2,3,4},g[]={.5,-2,1},v[3];
        velocity_update(u,d,g,v); eq(v[0],0);eq(v[1],8);eq(v[2],-1);
        // A pressure-relaxed gradient or another alpha_u would fail these values.
        FluxUpdate f{};f.density=2;f.old=10;f.alpha=.75;f.area[0]=3;
        f.velocity[0]=4;f.influence[0]=.5;f.compact_gradient[0]=5;f.reconstructed_gradient[0]=1;
        f.original_force[0]=3;f.reconstructed_force[0]=1;
        for(bool outlet:{false,true}) {
            eq(mass_flux_update(f,outlet,false),16); // fresh=24-12+6=18
            auto h=f; h.old=100;eq(mass_flux_update(h,outlet,false),38.5);
            h=f;h.alpha=1;eq(mass_flux_update(h,outlet,false),18);
            h=f;h.compact_gradient[0]+=.2;eq(mass_flux_update(h,outlet,false)-16,-.45);
            h=f;h.reconstructed_gradient[0]+=.2;eq(mass_flux_update(h,outlet,false)-16,.45);
            h=f;h.compact_gradient[0]=h.reconstructed_gradient[0];
            h.original_force[0]=h.reconstructed_force[0];eq(mass_flux_update(h,outlet,false),20.5);
        }
        eq(mass_flux_update(f,true,true),0);
        double prescribed[]={2,9,7},area[]={-3,0,0};
        eq(inlet_flux_update(2,prescribed,area,10,.75),-6.5);
        double p[]={1,3,5},shape[]={.5,.25,.25},a[]={0,3,4},mom[2];
        outlet_trace_moment(p,shape,a,false,mom);eq(mom[0],12.5);eq(mom[1],5);
        outlet_trace_moment(p,shape,a,true,mom);eq(mom[0],0);eq(mom[1],0);
        // Norms are summed per sample; opposed normals must not cancel the scalar area.
        double opposite[]={0,-3,-4},m1[2],m2[2];
        outlet_trace_moment(p,shape,a,false,m1);outlet_trace_moment(p,shape,opposite,false,m2);
        eq((m1[0]+m2[0])/(m1[1]+m2[1]),2.5);eq(m1[1]+m2[1],10);
        eq(outlet_trace_update(5,2,3,.05,99,false),3.9);
        eq(outlet_trace_update(5,2,3,1,99,false),2);
        eq(outlet_trace_update(5,2,3,.05,99,true),99);
        double face_area[]={1,0,0,1,0,0,1,0,0};
        double face_u[]={1,0,0,1,0,0,1,0,0},face_p[]={3,3,3},trace[]={2,2,2};
        const double variants[][3]={{2,-1,2},{-2,1,0},{0,0,0},{-1e-18,0,0}};
        for(int state:{0,1}) for(bool ignore:{false,true}) for(int k=0;k<4;++k) {
            int old_flags[]={state,state,state},flags[3];double result[3];
            outlet_reversal_update(variants[k],old_flags,face_u,face_p,trace,face_area,ignore,result,flags);
            for(int s=0;s<3;++s) {
                int expected=state;
                if(!ignore) expected=state ? 0 : (k==1 ? 1 : 0);
                eq(flags[s],expected);
                eq(result[s],!state && k==1 ? 0 : std::max(variants[k][s],0.0));
            }
        }
        // A wall reopens only when both the velocity and pressure tests permit it.
        int closed[]={1,1,1},flags[3];double result[3],positive[]={1,2,3};
        for(int cause=0;cause<2;++cause) {
            face_u[0]=face_u[3]=face_u[6]=cause==0 ? -1:1;
            trace[0]=trace[1]=trace[2]=cause==1 ? 4:2;
            outlet_reversal_update(positive,closed,face_u,face_p,trace,face_area,false,result,flags);
            for(int s=0;s<3;++s) {eq(result[s],0);eq(flags[s],1);}
        }
        std::cout<<"PASS: "<<checks<<" independent update algebra checks\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
