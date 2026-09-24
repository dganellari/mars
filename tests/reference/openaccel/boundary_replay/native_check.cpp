#include "mars_boundary_reference.hpp"
#include <algorithm>
#include <iostream>
#include <random>

int main() {
    using namespace mars::segregated;
    using Fn=void(*)(const BoundaryInput&,BoundaryOutput&);
    const Fn reference[]={boundary_reference::stage0,boundary_reference::stage1,boundary_reference::stage2,
                          boundary_reference::stage3,boundary_reference::stage4,boundary_reference::stage5};
    std::mt19937 gen(8731);std::uniform_real_distribution<double> value(-2,2);
    double worst[6]{};unsigned checks=0;
    for(int stage=0;stage<6;++stage)for(int trial=0;trial<200;++trial) {
        BoundaryInput x;x.stage=stage;
        const int count=stage==5?3:4, width=count*(stage<3?1:3);
        for(int s=0;s<3;++s) {
            x.face_nodes[s]=stage==5?s:(s+trial)%4;
            x.nearest[s]=x.face_nodes[s];x.opposing[s]=(trial+3)%4;
            x.reversal[s]=(trial+s)%3==0;
            x.density[s]=3+value(gen);x.viscosity[s]=3+value(gen);
            x.stored_flux[s]=value(gen);x.wall_coefficient[s]=3+value(gen);
            for(int j=0;j<3;++j) {
                x.area[3*s+j]=value(gen);x.shape[3*s+j]=s==j?11.0/18:7.0/36;
                x.boundary_velocity[3*s+j]=value(gen);
                x.influence_lhs[3*s+j]=3+value(gen);x.influence_rhs[3*s+j]=3+value(gen);
            }
            for(int n=0;n<4;++n)for(int j=0;j<3;++j)x.gradient[12*s+3*n+j]=value(gen);
        }
        for(int n=0;n<4;++n) {
            x.pressure[n]=value(gen);x.bc_multiplier[n]=n==(trial+3)%4?1:0;
            for(int j=0;j<3;++j){x.velocity[3*n+j]=value(gen);x.pressure_gradient[3*n+j]=value(gen);}
        }
        BoundaryOutput expected,actual;reference[stage](x,expected);boundary_block(x,actual);
        for(int group=0;group<2;++group) {
            const int length=group==0?width*width:width;
            const double* a=group==0?actual.lhs:actual.rhs;
            const double* b=group==0?expected.lhs:expected.rhs;
            double scale=1;for(int i=0;i<length;++i)scale=std::max(scale,std::abs(b[i]));
            for(int i=0;i<length;++i){++checks;double error=std::abs(a[i]-b[i])/scale;worst[stage]=std::max(worst[stage],error);
                if(!std::isfinite(a[i])||error>2e-13){std::cerr<<"FAIL: stage="<<stage<<" trial="<<trial<<" field="<<group<<" i="<<i<<" expected="<<b[i]<<" actual="<<a[i]<<'\n';return 1;}}
        }
    }
    for(int s=0;s<6;++s)std::cout<<"native stage="<<s<<" worst_scaled="<<worst[s]<<'\n';
    std::cout<<"PASS: "<<checks<<" checks against extracted reference arithmetic; not STK or GPU execution\n";
}
